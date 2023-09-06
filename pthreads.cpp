// Code written by Richard O. Lee and Christian Bienia
// Modified by Christian Fensch

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

#include "cha.hpp"
#include "topology.hpp"
#if defined(WIN32)
#define NOMINMAX
#include <windows.h>
#endif
#include <assert.h>
#include <float.h>
#include <math.h>
#include <pthread.h>

#include <map>
#include <mutex>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "cellpool.hpp"
#include "fluid.hpp"
#include "parsec_barrier.hpp"
#include "topology.hpp"

int getMostAccessedCHA(int tid1,
                       int tid2,
                       std::multiset<std::tuple<int, int, int, int>, std::greater<>> ranked_cha_access_count_per_pair,
                       Topology topo)
{
    int max = 0;
    std::vector<int> considered_chas;
    std::map<int, bool> considered_chas_flag;

    auto it = ranked_cha_access_count_per_pair.begin();
    while (it != ranked_cha_access_count_per_pair.end())
    {
        // std::pair<int, int> tid_pair(std::get<2>(*it), std::get<3>(*it));
        if ((std::get<2>(*it) == tid1 && std::get<3>(*it) == tid2) || (std::get<3>(*it) == tid1 && std::get<2>(*it) == tid2))
        {
            // SPDLOG_INFO("returning {}", std::get<1>(*it));
            max = std::get<0>(*it);
            considered_chas_flag[std::get<1>(*it)] = true;
            considered_chas.push_back(std::get<1>(*it));
            break;
        }
        it++;
    }

    // SPDLOG_INFO("communication between threads {} and {} uses the following chas the most", std::get<1>(*it), std::get<0>(*it), max);
    while (it != ranked_cha_access_count_per_pair.end())
    {
        // std::pair<int, int> tid_pair(std::get<2>(*it), std::get<3>(*it));
        if ((std::get<2>(*it) == tid1 && std::get<3>(*it) == tid2) || (std::get<3>(*it) == tid1 && std::get<2>(*it) == tid2))
        {
            // SPDLOG_INFO("returning {}", std::get<1>(*it));
            if (considered_chas_flag[std::get<1>(*it)] == false && std::get<0>(*it) > (0.9 * max))
            {
                // SPDLOG_INFO("cha {}, access count: {}, max: {}", std::get<1>(*it), std::get<0>(*it), max);
                considered_chas_flag[std::get<1>(*it)] = true;
                considered_chas.push_back(std::get<1>(*it));
            }
        }
        it++;
    }

    int x_total = 0;
    int y_total = 0;
    int cha_count = 0;
    // SPDLOG_INFO("communication between threads {} and {} involves the following CHAs:");
    for (auto it1 : considered_chas)
    {
        auto tile = topo.getTile(it1);
        x_total += tile.x;
        y_total += tile.y;
        cha_count++;
        // SPDLOG_INFO("cha {}, x: {}, y: {}", tile.cha, tile.x, tile.y);
    }

    int x_coord = x_total / cha_count;
    int y_coord = y_total / cha_count;
    auto tile = topo.getTile(x_coord, y_coord);
    // SPDLOG_INFO("the center of gravity is cha {}, x: {}, y: {}", tile.cha, tile.x, tile.y);
    // approximate the algorithm now
    return tile.cha;
}

void stick_this_thread_to_core(int core_id) {
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (core_id < 0 || core_id >= num_cores) {
        std::cerr << "error binding thread to core: " << core_id << '\n';
        // SPDLOG_ERROR("error binding thread to core {}!", core_id);
        return;
    }

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t current_thread = pthread_self();

    int res = pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    if (res == 0) {
        // std::cout << "thread bound to core " << core_id << std::endl;
        //        SPDLOG_INFO("Thread bound to core {} successfully.", core_id);
    } else {
        //        SPDLOG_ERROR("Error in binding this thread to core {}.", core_id);
    }
}

void assertRoot() {
    uid_t uid = getuid();
    if (uid == 0) {
        std::cout << "Running as root." << std::endl;
    } else {
        std::cerr << "Not running as root. Need root privileges to run the app. Exiting.\n";
        exit(EXIT_FAILURE);
    }
}

std::mutex map_mutex;
std::map<int, std::multiset<Cell *>>
    threadid_addresses_map;  // will use std::set_intersection. that is why I used a set. I picked multiset instead of
                             // set since if a thread pair communicates over same addresses multiple times, I want to
                             // take this into account.

#ifdef ENABLE_VISUALIZATION
#include "fluidview.hpp"
#endif

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

// Uncomment to add code to check that Courant–Friedrichs–Lewy condition is satisfied at runtime
//#define ENABLE_CFL_CHECK

////////////////////////////////////////////////////////////////////////////////

cellpool *pools;  // each thread has its private cell pool

fptype restParticlesPerMeter, h, hSq;
fptype densityCoeff, pressureCoeff, viscosityCoeff;

int nx, ny, nz;  // number of grid cells in each dimension
Vec3 delta;      // cell dimensions
int numParticles = 0;
int numCells = 0;
Cell *cells = 0;
Cell *cells2 = 0;
int *cnumPars = 0;
int *cnumPars2 = 0;
Cell **last_cells = NULL;  // helper array with pointers to last cell structure of "cells" array lists
#ifdef ENABLE_VISUALIZATION
Vec3 vMax(0.0, 0.0, 0.0);
Vec3 vMin(0.0, 0.0, 0.0);
#endif

int XDIVS = 1;  // number of partitions in X
int ZDIVS = 1;  // number of partitions in Z

#define NUM_GRIDS ((XDIVS) * (ZDIVS))
#define MUTEXES_PER_CELL 128
#define CELL_MUTEX_ID 0

struct Grid {
    union {
        struct {
            int sx, sy, sz;
            int ex, ey, ez;
        };
        unsigned char pp[CACHELINE_SIZE];
    };
} * grids;
bool *border;
pthread_attr_t attr;
pthread_t *thread;
pthread_mutex_t **mutex;    // used to lock cells in RebuildGrid and also particles in other functions
pthread_barrier_t barrier;  // global barrier used by all threads
#ifdef ENABLE_VISUALIZATION
pthread_barrier_t
    visualization_barrier;  // global barrier to separate (serial) visualization phase from (parallel) fluid simulation
#endif

typedef struct __thread_args {
    int tid;     // thread id, determines work partition
    int frames;  // number of frames to compute
} thread_args;   // arguments for threads

////////////////////////////////////////////////////////////////////////////////

/*
 * hmgweight
 *
 * Computes the hamming weight of x
 *
 * x      - input value
 * lsb    - if x!=0 position of smallest bit set, else -1
 *
 * return - the hamming weight
 */
unsigned int hmgweight(unsigned int x, int *lsb) {
    unsigned int weight = 0;
    unsigned int mask = 1;
    unsigned int count = 0;

    *lsb = -1;
    while (x > 0) {
        unsigned int temp;
        temp = (x & mask);
        if ((x & mask) == 1) {
            weight++;
            if (*lsb == -1) *lsb = count;
        }
        x >>= 1;
        count++;
    }

    return weight;
}

void InitSim(char const *fileName, unsigned int threadnum) {
    // Compute partitioning based on square root of number of threads
    // NOTE: Other partition sizes are possible as long as XDIVS * ZDIVS == threadnum,
    //       but communication is minimal (and hence optimal) if XDIVS == ZDIVS
    int lsb;
    if (hmgweight(threadnum, &lsb) != 1) {
        std::cerr << "Number of threads must be a power of 2" << std::endl;
        exit(1);
    }
    XDIVS = 1 << (lsb / 2);
    ZDIVS = 1 << (lsb / 2);
    if (XDIVS * ZDIVS != threadnum) XDIVS *= 2;
    assert(XDIVS * ZDIVS == threadnum);

    thread = new pthread_t[NUM_GRIDS];
    grids = new struct Grid[NUM_GRIDS];
    std::cout << "num grids: " << NUM_GRIDS << std::endl;
    assert(sizeof(Grid) <= CACHELINE_SIZE);  // as we put and aligh grid on the cacheline size to avoid false-sharing
                                             // if asserts fails - increase pp union member in Grid declarationi
                                             // and change this macro
    pools = new cellpool[NUM_GRIDS];

    // Load input particles
    std::cout << "Loading file \"" << fileName << "\"..." << std::endl;
    std::ifstream file(fileName, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file. Aborting." << std::endl;
        exit(1);
    }

    // Always use single precision float variables b/c file format uses single precision
    float restParticlesPerMeter_le;
    int numParticles_le;
    file.read((char *)&restParticlesPerMeter_le, FILE_SIZE_FLOAT);
    file.read((char *)&numParticles_le, FILE_SIZE_INT);
    if (!isLittleEndian()) {
        restParticlesPerMeter = bswap_float(restParticlesPerMeter_le);
        numParticles = bswap_int32(numParticles_le);
    } else {
        restParticlesPerMeter = restParticlesPerMeter_le;
        numParticles = numParticles_le;
    }
    for (int i = 0; i < NUM_GRIDS; i++) cellpool_init(&pools[i], numParticles / NUM_GRIDS);

    h = kernelRadiusMultiplier / restParticlesPerMeter;
    hSq = h * h;

#ifndef ENABLE_DOUBLE_PRECISION
    fptype coeff1 = 315.0 / (64.0 * pi * powf(h, 9.0));
    fptype coeff2 = 15.0 / (pi * powf(h, 6.0));
    fptype coeff3 = 45.0 / (pi * powf(h, 6.0));
#else
    fptype coeff1 = 315.0 / (64.0 * pi * pow(h, 9.0));
    fptype coeff2 = 15.0 / (pi * pow(h, 6.0));
    fptype coeff3 = 45.0 / (pi * pow(h, 6.0));
#endif  // ENABLE_DOUBLE_PRECISION
    fptype particleMass =
        0.5 * doubleRestDensity / (restParticlesPerMeter * restParticlesPerMeter * restParticlesPerMeter);
    densityCoeff = particleMass * coeff1;
    pressureCoeff = 3.0 * coeff2 * 0.50 * stiffnessPressure * particleMass;
    viscosityCoeff = viscosity * coeff3 * particleMass;

    Vec3 range = domainMax - domainMin;
    nx = (int)(range.x / h);
    ny = (int)(range.y / h);
    nz = (int)(range.z / h);
    assert(nx >= 1 && ny >= 1 && nz >= 1);
    numCells = nx * ny * nz;
    std::cout << "Number of cells: " << numCells << std::endl;
    delta.x = range.x / nx;
    delta.y = range.y / ny;
    delta.z = range.z / nz;
    assert(delta.x >= h && delta.y >= h && delta.z >= h);

    std::cout << "Grids steps over x, y, z: " << delta.x << " " << delta.y << " " << delta.z << std::endl;

    std::cout << "XDIVS: " << XDIVS << ", ZDIVS: " << ZDIVS << std::endl;

    assert(nx >= XDIVS && nz >= ZDIVS);
    int gi = 0;  // grid id.
    int sx, sz, ex, ez;
    ex = 0;
    for (int i = 0; i < XDIVS; ++i) {
        sx = ex;
        ex = (int)((fptype)(nx) / (fptype)(XDIVS) * (i + 1) + 0.5);
        assert(sx < ex);

        ez = 0;
        for (int j = 0; j < ZDIVS; ++j, ++gi) {
            sz = ez;
            ez = (int)((fptype)(nz) / (fptype)(ZDIVS) * (j + 1) + 0.5);
            assert(sz < ez);

            grids[gi].sx = sx;
            grids[gi].ex = ex;
            grids[gi].sy = 0;
            grids[gi].ey = ny;
            grids[gi].sz = sz;
            grids[gi].ez = ez;

            // std::cout << "gi: " << gi << ", sx: " << grids[gi].sx << ", ex: " << grids[gi].ex
            //           << ", sy: " << grids[gi].sy << ", ey: " << grids[gi].ey << ", sz: " << grids[gi].sz
            //           << ", ez: " << grids[gi].ez << std::endl;
        }
    }
    // std::exit(42);
    assert(gi == NUM_GRIDS);

    border = new bool[numCells];
    for (int i = 0; i < NUM_GRIDS; ++i)
        for (int iz = grids[i].sz; iz < grids[i].ez; ++iz)
            for (int iy = grids[i].sy; iy < grids[i].ey; ++iy)
                for (int ix = grids[i].sx; ix < grids[i].ex; ++ix) {
                    int index = (iz * ny + iy) * nx + ix;
                    border[index] = false;
                    for (int dk = -1; dk <= 1; ++dk) {
                        for (int dj = -1; dj <= 1; ++dj) {
                            for (int di = -1; di <= 1; ++di) {
                                int ci = ix + di;
                                int cj = iy + dj;
                                int ck = iz + dk;

                                if (ci < 0)
                                    ci = 0;
                                else if (ci > (nx - 1))
                                    ci = nx - 1;
                                if (cj < 0)
                                    cj = 0;
                                else if (cj > (ny - 1))
                                    cj = ny - 1;
                                if (ck < 0)
                                    ck = 0;
                                else if (ck > (nz - 1))
                                    ck = nz - 1;

                                if (ci < grids[i].sx || ci >= grids[i].ex || cj < grids[i].sy || cj >= grids[i].ey ||
                                    ck < grids[i].sz || ck >= grids[i].ez) {
                                    border[index] = true;
                                    break;
                                }
                            }  // for(int di = -1; di <= 1; ++di)
                            if (border[index]) break;
                        }  // for(int dj = -1; dj <= 1; ++dj)
                        if (border[index]) break;
                    }  // for(int dk = -1; dk <= 1; ++dk)
                }

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    mutex = new pthread_mutex_t *[numCells];
    for (int i = 0; i < numCells; ++i) {
        assert(CELL_MUTEX_ID < MUTEXES_PER_CELL);
        int n = (border[i] ? MUTEXES_PER_CELL : CELL_MUTEX_ID + 1);
        mutex[i] = new pthread_mutex_t[n];
        for (int j = 0; j < n; ++j) pthread_mutex_init(&mutex[i][j], NULL);
    }
    pthread_barrier_init(&barrier, NULL, NUM_GRIDS);
#ifdef ENABLE_VISUALIZATION
    // visualization barrier is used by all NUM_GRIDS worker threads and 1 master thread
    pthread_barrier_init(&visualization_barrier, NULL, NUM_GRIDS + 1);
#endif
    // make sure Cell structure is multiple of estiamted cache line size
    assert(sizeof(Cell) % CACHELINE_SIZE == 0);
    // make sure helper Cell structure is in sync with real Cell structure
    assert(offsetof(struct Cell_aux, padding) == offsetof(struct Cell, padding));

#if defined(WIN32)
    cells = (struct Cell *)_aligned_malloc(sizeof(struct Cell) * numCells, CACHELINE_SIZE);
    cells2 = (struct Cell *)_aligned_malloc(sizeof(struct Cell) * numCells, CACHELINE_SIZE);
    cnumPars = (int *)_aligned_malloc(sizeof(int) * numCells, CACHELINE_SIZE);
    cnumPars2 = (int *)_aligned_malloc(sizeof(int) * numCells, CACHELINE_SIZE);
    last_cells = (struct Cell **)_aligned_malloc(sizeof(struct Cell *) * numCells, CACHELINE_SIZE);
    assert((cells != NULL) && (cells2 != NULL) && (cnumPars != NULL) && (cnumPars2 != NULL) && (last_cells != NULL));
#elif defined(SPARC_SOLARIS)
    cells = (Cell *)memalign(CACHELINE_SIZE, sizeof(struct Cell) * numCells);
    cells2 = (Cell *)memalign(CACHELINE_SIZE, sizeof(struct Cell) * numCells);
    cnumPars = (int *)memalign(CACHELINE_SIZE, sizeof(int) * numCells);
    cnumPars2 = (int *)memalign(CACHELINE_SIZE, sizeof(int) * numCells);
    last_cells = (Cell **)memalign(CACHELINE_SIZE, sizeof(struct Cell *) * numCells);
    assert((cells != 0) && (cells2 != 0) && (cnumPars != 0) && (cnumPars2 != 0) && (last_cells != 0));
#else
    int rv0 = posix_memalign((void **)(&cells), CACHELINE_SIZE, sizeof(struct Cell) * numCells);
    int rv1 = posix_memalign((void **)(&cells2), CACHELINE_SIZE, sizeof(struct Cell) * numCells);
    int rv2 = posix_memalign((void **)(&cnumPars), CACHELINE_SIZE, sizeof(int) * numCells);
    int rv3 = posix_memalign((void **)(&cnumPars2), CACHELINE_SIZE, sizeof(int) * numCells);
    int rv4 = posix_memalign((void **)(&last_cells), CACHELINE_SIZE, sizeof(struct Cell *) * numCells);
    assert((rv0 == 0) && (rv1 == 0) && (rv2 == 0) && (rv3 == 0) && (rv4 == 0));
#endif

    // because cells and cells2 are not allocated via new
    // we construct them here
    for (int i = 0; i < numCells; ++i) {
        new (&cells[i]) Cell;
        new (&cells2[i]) Cell;
    }

    memset(cnumPars, 0, numCells * sizeof(int));

    // Always use single precision float variables b/c file format uses single precision float
    int pool_id = 0;
    float px, py, pz, hvx, hvy, hvz, vx, vy, vz;
    for (int i = 0; i < numParticles; ++i) {
        file.read((char *)&px, FILE_SIZE_FLOAT);
        file.read((char *)&py, FILE_SIZE_FLOAT);
        file.read((char *)&pz, FILE_SIZE_FLOAT);
        file.read((char *)&hvx, FILE_SIZE_FLOAT);
        file.read((char *)&hvy, FILE_SIZE_FLOAT);
        file.read((char *)&hvz, FILE_SIZE_FLOAT);
        file.read((char *)&vx, FILE_SIZE_FLOAT);
        file.read((char *)&vy, FILE_SIZE_FLOAT);
        file.read((char *)&vz, FILE_SIZE_FLOAT);
        if (!isLittleEndian()) {
            px = bswap_float(px);
            py = bswap_float(py);
            pz = bswap_float(pz);
            hvx = bswap_float(hvx);
            hvy = bswap_float(hvy);
            hvz = bswap_float(hvz);
            vx = bswap_float(vx);
            vy = bswap_float(vy);
            vz = bswap_float(vz);
        }

        int ci = (int)((px - domainMin.x) / delta.x);
        int cj = (int)((py - domainMin.y) / delta.y);
        int ck = (int)((pz - domainMin.z) / delta.z);

        if (ci < 0)
            ci = 0;
        else if (ci > (nx - 1))
            ci = nx - 1;
        if (cj < 0)
            cj = 0;
        else if (cj > (ny - 1))
            cj = ny - 1;
        if (ck < 0)
            ck = 0;
        else if (ck > (nz - 1))
            ck = nz - 1;

        int index = (ck * ny + cj) * nx + ci;
        Cell *cell = &cells[index];

        // go to last cell structure in list
        int np = cnumPars[index];
        while (np > PARTICLES_PER_CELL) {
            cell = cell->next;
            np = np - PARTICLES_PER_CELL;
        }
        // add another cell structure if everything full
        if ((np % PARTICLES_PER_CELL == 0) && (cnumPars[index] != 0)) {
            // Get cells from pools in round-robin fashion to balance load during parallel phase
            cell->next = cellpool_getcell(&pools[pool_id]);
            pool_id = (pool_id + 1) % NUM_GRIDS;
            cell = cell->next;
            np = np - PARTICLES_PER_CELL;
        }

        cell->p[np].x = px;
        cell->p[np].y = py;
        cell->p[np].z = pz;
        cell->hv[np].x = hvx;
        cell->hv[np].y = hvy;
        cell->hv[np].z = hvz;
        cell->v[np].x = vx;
        cell->v[np].y = vy;
        cell->v[np].z = vz;
#ifdef ENABLE_VISUALIZATION
        vMin.x = std::min(vMin.x, cell->v[np].x);
        vMax.x = std::max(vMax.x, cell->v[np].x);
        vMin.y = std::min(vMin.y, cell->v[np].y);
        vMax.y = std::max(vMax.y, cell->v[np].y);
        vMin.z = std::min(vMin.z, cell->v[np].z);
        vMax.z = std::max(vMax.z, cell->v[np].z);
#endif
        ++cnumPars[index];
    }

    std::cout << "Number of particles: " << numParticles << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void SaveFile(char const *fileName) {
    std::cout << "Saving file \"" << fileName << "\"..." << std::endl;

    std::ofstream file(fileName, std::ios::binary);
    assert(file);

    // Always use single precision float variables b/c file format uses single precision
    if (!isLittleEndian()) {
        float restParticlesPerMeter_le;
        int numParticles_le;

        restParticlesPerMeter_le = bswap_float((float)restParticlesPerMeter);
        numParticles_le = bswap_int32(numParticles);
        file.write((char *)&restParticlesPerMeter_le, FILE_SIZE_FLOAT);
        file.write((char *)&numParticles_le, FILE_SIZE_INT);
    } else {
        file.write((char *)&restParticlesPerMeter, FILE_SIZE_FLOAT);
        file.write((char *)&numParticles, FILE_SIZE_INT);
    }

    int count = 0;
    for (int i = 0; i < numCells; ++i) {
        Cell *cell = &cells[i];
        int np = cnumPars[i];
        for (int j = 0; j < np; ++j) {
            // Always use single precision float variables b/c file format uses single precision
            float px, py, pz, hvx, hvy, hvz, vx, vy, vz;
            if (!isLittleEndian()) {
                px = bswap_float((float)(cell->p[j % PARTICLES_PER_CELL].x));
                py = bswap_float((float)(cell->p[j % PARTICLES_PER_CELL].y));
                pz = bswap_float((float)(cell->p[j % PARTICLES_PER_CELL].z));
                hvx = bswap_float((float)(cell->hv[j % PARTICLES_PER_CELL].x));
                hvy = bswap_float((float)(cell->hv[j % PARTICLES_PER_CELL].y));
                hvz = bswap_float((float)(cell->hv[j % PARTICLES_PER_CELL].z));
                vx = bswap_float((float)(cell->v[j % PARTICLES_PER_CELL].x));
                vy = bswap_float((float)(cell->v[j % PARTICLES_PER_CELL].y));
                vz = bswap_float((float)(cell->v[j % PARTICLES_PER_CELL].z));
            } else {
                px = (float)(cell->p[j % PARTICLES_PER_CELL].x);
                py = (float)(cell->p[j % PARTICLES_PER_CELL].y);
                pz = (float)(cell->p[j % PARTICLES_PER_CELL].z);
                hvx = (float)(cell->hv[j % PARTICLES_PER_CELL].x);
                hvy = (float)(cell->hv[j % PARTICLES_PER_CELL].y);
                hvz = (float)(cell->hv[j % PARTICLES_PER_CELL].z);
                vx = (float)(cell->v[j % PARTICLES_PER_CELL].x);
                vy = (float)(cell->v[j % PARTICLES_PER_CELL].y);
                vz = (float)(cell->v[j % PARTICLES_PER_CELL].z);
            }
            file.write((char *)&px, FILE_SIZE_FLOAT);
            file.write((char *)&py, FILE_SIZE_FLOAT);
            file.write((char *)&pz, FILE_SIZE_FLOAT);
            file.write((char *)&hvx, FILE_SIZE_FLOAT);
            file.write((char *)&hvy, FILE_SIZE_FLOAT);
            file.write((char *)&hvz, FILE_SIZE_FLOAT);
            file.write((char *)&vx, FILE_SIZE_FLOAT);
            file.write((char *)&vy, FILE_SIZE_FLOAT);
            file.write((char *)&vz, FILE_SIZE_FLOAT);
            ++count;

            // move pointer to next cell in list if end of array is reached
            if (j % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                cell = cell->next;
            }
        }
    }
    assert(count == numParticles);
}

////////////////////////////////////////////////////////////////////////////////

void CleanUpSim() {
    // first return extended cells to cell pools
    for (int i = 0; i < numCells; ++i) {
        Cell &cell = cells[i];
        while (cell.next) {
            Cell *temp = cell.next;
            cell.next = temp->next;
            cellpool_returncell(&pools[0], temp);
        }
    }
    // now return cell pools
    // NOTE: Cells from cell pools can migrate to different pools during the parallel phase.
    //      This is no problem as long as all cell pools are destroyed together. Each pool
    //      uses its internal meta information to free exactly the cells which it allocated
    //      itself. This guarantees that all allocated cells will be freed but it might
    //      render other cell pools unusable so they also have to be destroyed.
    for (int i = 0; i < NUM_GRIDS; i++) cellpool_destroy(&pools[i]);
    pthread_attr_destroy(&attr);

    for (int i = 0; i < numCells; ++i) {
        assert(CELL_MUTEX_ID < MUTEXES_PER_CELL);
        int n = (border[i] ? MUTEXES_PER_CELL : CELL_MUTEX_ID + 1);
        for (int j = 0; j < n; ++j) pthread_mutex_destroy(&mutex[i][j]);
        delete[] mutex[i];
    }
    delete[] mutex;
    pthread_barrier_destroy(&barrier);
#ifdef ENABLE_VISUALIZATION
    pthread_barrier_destroy(&visualization_barrier);
#endif

    delete[] border;

#if defined(WIN32)
    _aligned_free(cells);
    _aligned_free(cells2);
    _aligned_free(cnumPars);
    _aligned_free(cnumPars2);
    _aligned_free(last_cells);
#else
    free(cells);
    free(cells2);
    free(cnumPars);
    free(cnumPars2);
    free(last_cells);
#endif
    delete[] thread;
    delete[] grids;
}

////////////////////////////////////////////////////////////////////////////////
void ClearParticlesMT(int tid) {
    // std::cout << __PRETTY_FUNCTION__ << std::endl;
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz) {
        // std::cout << "tid: " << tid << ", sz: " << grids[tid].sz << ", ez: " << grids[tid].ez << std::endl;
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy) {
            // std::cout << "tid: " << tid << ", sy: " << grids[tid].sy << ", ey: " << grids[tid].ey << std::endl;
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                // std::cout << "tid: " << tid << ", sx: " << grids[tid].sx << ", ex: " << grids[tid].ex << std::endl;
                int index = (iz * ny + iy) * nx + ix;
                cnumPars[index] = 0;
                cells[index].next = NULL;
                last_cells[index] = &cells[index];
            }
        }
    }
    // std::exit(42);
}

////////////////////////////////////////////////////////////////////////////////

void RebuildGridMT(int tid) {
    // Note, in parallel versions the below swaps
    // occure outside RebuildGrid()
    // swap src and dest arrays with particles
    //   std::swap(cells, cells2);
    // swap src and dest arrays with counts of particles
    //  std::swap(cnumPars, cnumPars2);

    // iterate through source cell lists
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                int index2 = (iz * ny + iy) * nx + ix;

                Cell *cell2 = &cells2[index2];

                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell2);
                }

                int np2 = cnumPars2[index2];
                // iterate through source particles
                for (int j = 0; j < np2; ++j) {
                    // get destination for source particle
                    int ci = (int)((cell2->p[j % PARTICLES_PER_CELL].x - domainMin.x) / delta.x);
                    int cj = (int)((cell2->p[j % PARTICLES_PER_CELL].y - domainMin.y) / delta.y);
                    int ck = (int)((cell2->p[j % PARTICLES_PER_CELL].z - domainMin.z) / delta.z);

                    if (ci < 0)
                        ci = 0;
                    else if (ci > (nx - 1))
                        ci = nx - 1;
                    if (cj < 0)
                        cj = 0;
                    else if (cj > (ny - 1))
                        cj = ny - 1;
                    if (ck < 0)
                        ck = 0;
                    else if (ck > (nz - 1))
                        ck = nz - 1;
#if 0
		  assert(ci>=ix-1);
		  assert(ci<=ix+1);
		  assert(cj>=iy-1);
		  assert(cj<=iy+1);
		  assert(ck>=iz-1);
		  assert(ck<=iz+1);
#endif
#ifdef ENABLE_CFL_CHECK
                    // check that source cell is a neighbor of destination cell
                    bool cfl_cond_satisfied = false;
                    for (int di = -1; di <= 1; ++di)
                        for (int dj = -1; dj <= 1; ++dj)
                            for (int dk = -1; dk <= 1; ++dk) {
                                int ii = ci + di;
                                int jj = cj + dj;
                                int kk = ck + dk;
                                if (ii >= 0 && ii < nx && jj >= 0 && jj < ny && kk >= 0 && kk < nz) {
                                    int index = (kk * ny + jj) * nx + ii;
                                    if (index == index2) {
                                        cfl_cond_satisfied = true;
                                        break;
                                    }
                                }
                            }
                    if (!cfl_cond_satisfied) {
                        std::cerr << "FATAL ERROR: Courant–Friedrichs–Lewy condition not satisfied." << std::endl;
                        exit(1);
                    }
#endif  // ENABLE_CFL_CHECK

                    int index = (ck * ny + cj) * nx + ci;
                    // this assumes that particles cannot travel more than one grid cell per time step
                    if (border[index]) pthread_mutex_lock(&mutex[index][CELL_MUTEX_ID]);
                    Cell *cell = last_cells[index];
                    {
                        std::lock_guard lock(map_mutex);
                        threadid_addresses_map[tid].insert(cell);
                    }
                    int np = cnumPars[index];

                    // add another cell structure if everything full
                    if ((np % PARTICLES_PER_CELL == 0) && (cnumPars[index] != 0)) {
                        cell->next = cellpool_getcell(&pools[tid]);
                        cell = cell->next;
                        last_cells[index] = cell;
                    }
                    ++cnumPars[index];
                    if (border[index]) pthread_mutex_unlock(&mutex[index][CELL_MUTEX_ID]);

                    // copy source to destination particle

                    cell->p[np % PARTICLES_PER_CELL] = cell2->p[j % PARTICLES_PER_CELL];
                    cell->hv[np % PARTICLES_PER_CELL] = cell2->hv[j % PARTICLES_PER_CELL];
                    cell->v[np % PARTICLES_PER_CELL] = cell2->v[j % PARTICLES_PER_CELL];
                    // move pointer to next source cell in list if end of array is reached
                    if (j % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                        Cell *temp = cell2;

                        // {
                        //     std::lock_guard lock(map_mutex);
                        //     address_thread_ids_map[temp].insert(tid);
                        // }

                        cell2 = cell2->next;
                        // return cells to pool that are not statically allocated head of lists
                        if (temp != &cells2[index2]) {
                            // NOTE: This is thread-safe because temp and pool are thread-private, no need to
                            // synchronize
                            cellpool_returncell(&pools[tid], temp);
                        }
                    }
                }  // for(int j = 0; j < np2; ++j)
                // return cells to pool that are not statically allocated head of lists
                if ((cell2 != NULL) && (cell2 != &cells2[index2])) {
                    cellpool_returncell(&pools[tid], cell2);
                }
            }
}

////////////////////////////////////////////////////////////////////////////////

int InitNeighCellList(int ci, int cj, int ck, int *neighCells) {
    int numNeighCells = 0;

    // have the nearest particles first -> help branch prediction
    int my_index = (ck * ny + cj) * nx + ci;
    neighCells[numNeighCells] = my_index;
    ++numNeighCells;

    for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
            for (int dk = -1; dk <= 1; ++dk) {
                int ii = ci + di;
                int jj = cj + dj;
                int kk = ck + dk;
                if (ii >= 0 && ii < nx && jj >= 0 && jj < ny && kk >= 0 && kk < nz) {
                    int index = (kk * ny + jj) * nx + ii;
                    if ((index < my_index) && (cnumPars[index] != 0)) {
                        neighCells[numNeighCells] = index;
                        ++numNeighCells;
                    }
                }
            }
    return numNeighCells;
}

////////////////////////////////////////////////////////////////////////////////

void InitDensitiesAndForcesMT(int tid) {
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                int index = (iz * ny + iy) * nx + ix;
                Cell *cell = &cells[index];

                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }

                int np = cnumPars[index];
                for (int j = 0; j < np; ++j) {
                    cell->density[j % PARTICLES_PER_CELL] = 0.0;
                    cell->a[j % PARTICLES_PER_CELL] = externalAcceleration;
                    // move pointer to next cell in list if end of array is reached
                    if (j % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
}

////////////////////////////////////////////////////////////////////////////////

void ComputeDensitiesMT(int tid) {
    int neighCells[3 * 3 * 3];

    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                int index = (iz * ny + iy) * nx + ix;
                int np = cnumPars[index];
                if (np == 0) continue;

                int numNeighCells = InitNeighCellList(ix, iy, iz, neighCells);

                Cell *cell = &cells[index];
                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }
                for (int ipar = 0; ipar < np; ++ipar) {
                    for (int inc = 0; inc < numNeighCells; ++inc) {
                        int indexNeigh = neighCells[inc];
                        Cell *neigh = &cells[indexNeigh];
                        {
                            std::lock_guard lock(map_mutex);
                            threadid_addresses_map[tid].insert(neigh);
                        }
                        int numNeighPars = cnumPars[indexNeigh];
                        for (int iparNeigh = 0; iparNeigh < numNeighPars; ++iparNeigh) {
                            // Check address to make sure densities are computed only once per pair
                            if (&neigh->p[iparNeigh % PARTICLES_PER_CELL] < &cell->p[ipar % PARTICLES_PER_CELL]) {
                                fptype distSq =
                                    (cell->p[ipar % PARTICLES_PER_CELL] - neigh->p[iparNeigh % PARTICLES_PER_CELL])
                                        .GetLengthSq();
                                if (distSq < hSq) {
                                    fptype t = hSq - distSq;
                                    fptype tc = t * t * t;

                                    if (border[index]) {
                                        pthread_mutex_lock(&mutex[index][ipar % MUTEXES_PER_CELL]);
                                        cell->density[ipar % PARTICLES_PER_CELL] += tc;
                                        pthread_mutex_unlock(&mutex[index][ipar % MUTEXES_PER_CELL]);
                                    } else
                                        cell->density[ipar % PARTICLES_PER_CELL] += tc;

                                    if (border[indexNeigh]) {
                                        pthread_mutex_lock(&mutex[indexNeigh][iparNeigh % MUTEXES_PER_CELL]);
                                        neigh->density[iparNeigh % PARTICLES_PER_CELL] += tc;
                                        pthread_mutex_unlock(&mutex[indexNeigh][iparNeigh % MUTEXES_PER_CELL]);
                                    } else
                                        neigh->density[iparNeigh % PARTICLES_PER_CELL] += tc;
                                }
                            }
                            // move pointer to next cell in list if end of array is reached
                            if (iparNeigh % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                                neigh = neigh->next;
                            }
                        }
                    }
                    // move pointer to next cell in list if end of array is reached
                    if (ipar % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
}

////////////////////////////////////////////////////////////////////////////////

void ComputeDensities2MT(int tid) {
    const fptype tc = hSq * hSq * hSq;
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                int index = (iz * ny + iy) * nx + ix;
                Cell *cell = &cells[index];
                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }
                int np = cnumPars[index];
                for (int j = 0; j < np; ++j) {
                    cell->density[j % PARTICLES_PER_CELL] += tc;
                    cell->density[j % PARTICLES_PER_CELL] *= densityCoeff;
                    // move pointer to next cell in list if end of array is reached
                    if (j % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
}

////////////////////////////////////////////////////////////////////////////////

void ComputeForcesMT(int tid) {
    int neighCells[3 * 3 * 3];

    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                int index = (iz * ny + iy) * nx + ix;
                int np = cnumPars[index];
                if (np == 0) continue;

                int numNeighCells = InitNeighCellList(ix, iy, iz, neighCells);

                Cell *cell = &cells[index];
                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }
                for (int ipar = 0; ipar < np; ++ipar) {
                    for (int inc = 0; inc < numNeighCells; ++inc) {
                        int indexNeigh = neighCells[inc];
                        Cell *neigh = &cells[indexNeigh];
                        {
                            std::lock_guard lock(map_mutex);
                            threadid_addresses_map[tid].insert(neigh);
                        }
                        int numNeighPars = cnumPars[indexNeigh];
                        for (int iparNeigh = 0; iparNeigh < numNeighPars; ++iparNeigh) {
                            // Check address to make sure forces are computed only once per pair
                            if (&neigh->p[iparNeigh % PARTICLES_PER_CELL] < &cell->p[ipar % PARTICLES_PER_CELL]) {
                                Vec3 disp =
                                    cell->p[ipar % PARTICLES_PER_CELL] - neigh->p[iparNeigh % PARTICLES_PER_CELL];
                                fptype distSq = disp.GetLengthSq();
                                if (distSq < hSq) {
#ifndef ENABLE_DOUBLE_PRECISION
                                    fptype dist = sqrtf(std::max(distSq, (fptype)1e-12));
#else
                                    fptype dist = sqrt(std::max(distSq, 1e-12));
#endif  // ENABLE_DOUBLE_PRECISION
                                    fptype hmr = h - dist;

                                    Vec3 acc = disp * pressureCoeff * (hmr * hmr / dist) *
                                               (cell->density[ipar % PARTICLES_PER_CELL] +
                                                neigh->density[iparNeigh % PARTICLES_PER_CELL] - doubleRestDensity);
                                    acc += (neigh->v[iparNeigh % PARTICLES_PER_CELL] -
                                            cell->v[ipar % PARTICLES_PER_CELL]) *
                                           viscosityCoeff * hmr;
                                    acc /= cell->density[ipar % PARTICLES_PER_CELL] *
                                           neigh->density[iparNeigh % PARTICLES_PER_CELL];

                                    if (border[index]) {
                                        pthread_mutex_lock(&mutex[index][ipar % MUTEXES_PER_CELL]);
                                        cell->a[ipar % PARTICLES_PER_CELL] += acc;
                                        pthread_mutex_unlock(&mutex[index][ipar % MUTEXES_PER_CELL]);
                                    } else
                                        cell->a[ipar % PARTICLES_PER_CELL] += acc;

                                    if (border[indexNeigh]) {
                                        pthread_mutex_lock(&mutex[indexNeigh][iparNeigh % MUTEXES_PER_CELL]);
                                        neigh->a[iparNeigh % PARTICLES_PER_CELL] -= acc;
                                        pthread_mutex_unlock(&mutex[indexNeigh][iparNeigh % MUTEXES_PER_CELL]);
                                    } else
                                        neigh->a[iparNeigh % PARTICLES_PER_CELL] -= acc;
                                }
                            }
                            // move pointer to next cell in list if end of array is reached
                            if (iparNeigh % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                                neigh = neigh->next;
                            }
                        }
                    }
                    // move pointer to next cell in list if end of array is reached
                    if (ipar % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
}

////////////////////////////////////////////////////////////////////////////////

// ProcessCollisions() with container walls
// Under the assumptions that
// a) a particle will not penetrate a wall
// b) a particle will not migrate further than once cell
// c) the parSize is smaller than a cell
// then only the particles at the perimiters may be influenced by the walls
#if 0
void ProcessCollisionsMT(int tid)
{
  for(int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
    for(int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
      for(int ix = grids[tid].sx; ix < grids[tid].ex; ++ix)
      {
        int index = (iz*ny + iy)*nx + ix;
        Cell *cell = &cells[index];
        int np = cnumPars[index];
        for(int j = 0; j < np; ++j)
        {
          Vec3 pos = cell->p[j % PARTICLES_PER_CELL] + cell->hv[j % PARTICLES_PER_CELL] * timeStep;

          fptype diff = parSize - (pos.x - domainMin.x);
          if(diff > epsilon)
            cell->a[j % PARTICLES_PER_CELL].x += stiffnessCollisions*diff - damping*cell->v[j % PARTICLES_PER_CELL].x;

          diff = parSize - (domainMax.x - pos.x);
          if(diff > epsilon)
            cell->a[j % PARTICLES_PER_CELL].x -= stiffnessCollisions*diff + damping*cell->v[j % PARTICLES_PER_CELL].x;

          diff = parSize - (pos.y - domainMin.y);
          if(diff > epsilon)
            cell->a[j % PARTICLES_PER_CELL].y += stiffnessCollisions*diff - damping*cell->v[j % PARTICLES_PER_CELL].y;

          diff = parSize - (domainMax.y - pos.y);
          if(diff > epsilon)
            cell->a[j % PARTICLES_PER_CELL].y -= stiffnessCollisions*diff + damping*cell->v[j % PARTICLES_PER_CELL].y;

          diff = parSize - (pos.z - domainMin.z);
          if(diff > epsilon)
            cell->a[j % PARTICLES_PER_CELL].z += stiffnessCollisions*diff - damping*cell->v[j % PARTICLES_PER_CELL].z;

          diff = parSize - (domainMax.z - pos.z);
          if(diff > epsilon)
            cell->a[j % PARTICLES_PER_CELL].z -= stiffnessCollisions*diff + damping*cell->v[j % PARTICLES_PER_CELL].z;

          //move pointer to next cell in list if end of array is reached
          if(j % PARTICLES_PER_CELL == PARTICLES_PER_CELL-1) {
            cell = cell->next;
          }
        }
      }
}
#else
void ProcessCollisionsMT(int tid) {
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz) {
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy) {
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                if (!((ix == 0) || (iy == 0) || (iz == 0) || (ix == (nx - 1)) || (iy == (ny - 1)) == (iz == (nz - 1))))
                    continue;  // not on domain wall
                int index = (iz * ny + iy) * nx + ix;
                Cell *cell = &cells[index];
                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }
                int np = cnumPars[index];
                for (int j = 0; j < np; ++j) {
                    int ji = j % PARTICLES_PER_CELL;
                    Vec3 pos = cell->p[ji] + cell->hv[ji] * timeStep;

                    if (ix == 0) {
                        fptype diff = parSize - (pos.x - domainMin.x);
                        if (diff > epsilon) cell->a[ji].x += stiffnessCollisions * diff - damping * cell->v[ji].x;
                    }
                    if (ix == (nx - 1)) {
                        fptype diff = parSize - (domainMax.x - pos.x);
                        if (diff > epsilon) cell->a[ji].x -= stiffnessCollisions * diff + damping * cell->v[ji].x;
                    }
                    if (iy == 0) {
                        fptype diff = parSize - (pos.y - domainMin.y);
                        if (diff > epsilon) cell->a[ji].y += stiffnessCollisions * diff - damping * cell->v[ji].y;
                    }
                    if (iy == (ny - 1)) {
                        fptype diff = parSize - (domainMax.y - pos.y);
                        if (diff > epsilon) cell->a[ji].y -= stiffnessCollisions * diff + damping * cell->v[ji].y;
                    }
                    if (iz == 0) {
                        fptype diff = parSize - (pos.z - domainMin.z);
                        if (diff > epsilon) cell->a[ji].z += stiffnessCollisions * diff - damping * cell->v[ji].z;
                    }
                    if (iz == (nz - 1)) {
                        fptype diff = parSize - (domainMax.z - pos.z);
                        if (diff > epsilon) cell->a[ji].z -= stiffnessCollisions * diff + damping * cell->v[ji].z;
                    }
                    // move pointer to next cell in list if end of array is reached
                    if (ji == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
        }
    }
}
#endif

#define USE_ImpeneratableWall
#if defined(USE_ImpeneratableWall)
void ProcessCollisions2MT(int tid) {
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz) {
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy) {
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
#if 0
// Chris, the following test should be valid
// *** provided that a particle does not migrate more than 1 cell
// *** per integration step. This does not appear to be the case
// *** in the pthreads version. Serial version it seems to be OK
	    if(!((ix==0)||(iy==0)||(iz==0)||(ix==(nx-1))||(iy==(ny-1))==(iz==(nz-1))))
			continue;	// not on domain wall
#endif
                int index = (iz * ny + iy) * nx + ix;
                Cell *cell = &cells[index];
                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }
                int np = cnumPars[index];
                for (int j = 0; j < np; ++j) {
                    int ji = j % PARTICLES_PER_CELL;
                    Vec3 pos = cell->p[ji];

                    if (ix == 0) {
                        fptype diff = pos.x - domainMin.x;
                        if (diff < Zero) {
                            cell->p[ji].x = domainMin.x - diff;
                            cell->v[ji].x = -cell->v[ji].x;
                            cell->hv[ji].x = -cell->hv[ji].x;
                        }
                    }
                    if (ix == (nx - 1)) {
                        fptype diff = domainMax.x - pos.x;
                        if (diff < Zero) {
                            cell->p[ji].x = domainMax.x + diff;
                            cell->v[ji].x = -cell->v[ji].x;
                            cell->hv[ji].x = -cell->hv[ji].x;
                        }
                    }
                    if (iy == 0) {
                        fptype diff = pos.y - domainMin.y;
                        if (diff < Zero) {
                            cell->p[ji].y = domainMin.y - diff;
                            cell->v[ji].y = -cell->v[ji].y;
                            cell->hv[ji].y = -cell->hv[ji].y;
                        }
                    }
                    if (iy == (ny - 1)) {
                        fptype diff = domainMax.y - pos.y;
                        if (diff < Zero) {
                            cell->p[ji].y = domainMax.y + diff;
                            cell->v[ji].y = -cell->v[ji].y;
                            cell->hv[ji].y = -cell->hv[ji].y;
                        }
                    }
                    if (iz == 0) {
                        fptype diff = pos.z - domainMin.z;
                        if (diff < Zero) {
                            cell->p[ji].z = domainMin.z - diff;
                            cell->v[ji].z = -cell->v[ji].z;
                            cell->hv[ji].z = -cell->hv[ji].z;
                        }
                    }
                    if (iz == (nz - 1)) {
                        fptype diff = domainMax.z - pos.z;
                        if (diff < Zero) {
                            cell->p[ji].z = domainMax.z + diff;
                            cell->v[ji].z = -cell->v[ji].z;
                            cell->hv[ji].z = -cell->hv[ji].z;
                        }
                    }
                    // move pointer to next cell in list if end of array is reached
                    if (ji == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
        }
    }
}
#endif

////////////////////////////////////////////////////////////////////////////////

void AdvanceParticlesMT(int tid) {
    for (int iz = grids[tid].sz; iz < grids[tid].ez; ++iz)
        for (int iy = grids[tid].sy; iy < grids[tid].ey; ++iy)
            for (int ix = grids[tid].sx; ix < grids[tid].ex; ++ix) {
                int index = (iz * ny + iy) * nx + ix;
                Cell *cell = &cells[index];
                {
                    std::lock_guard lock(map_mutex);
                    threadid_addresses_map[tid].insert(cell);
                }
                int np = cnumPars[index];
                for (int j = 0; j < np; ++j) {
                    Vec3 v_half = cell->hv[j % PARTICLES_PER_CELL] + cell->a[j % PARTICLES_PER_CELL] * timeStep;
#if defined(USE_ImpeneratableWall)
                    // N.B. The integration of the position can place the particle
                    // outside the domain. Although we could place a test in this loop
                    // we would be unnecessarily testing particles on interior cells.
                    // Therefore, to reduce the amount of computations we make a later
                    // pass on the perimiter cells to account for particle migration
                    // beyond domain
#endif
                    cell->p[j % PARTICLES_PER_CELL] += v_half * timeStep;
                    cell->v[j % PARTICLES_PER_CELL] = cell->hv[j % PARTICLES_PER_CELL] + v_half;
                    cell->v[j % PARTICLES_PER_CELL] *= 0.5;
                    cell->hv[j % PARTICLES_PER_CELL] = v_half;

                    // move pointer to next cell in list if end of array is reached
                    if (j % PARTICLES_PER_CELL == PARTICLES_PER_CELL - 1) {
                        cell = cell->next;
                    }
                }
            }
}

////////////////////////////////////////////////////////////////////////////////

void AdvanceFrameMT(int tid) {
    // std::cout << __PRETTY_FUNCTION__ << std::endl;
    // swap src and dest arrays with particles
    if (tid == 0) {
        std::swap(cells, cells2);
        std::swap(cnumPars, cnumPars2);
    }
    pthread_barrier_wait(&barrier);

    ClearParticlesMT(tid);
    pthread_barrier_wait(&barrier);
    RebuildGridMT(tid);
    pthread_barrier_wait(&barrier);
    InitDensitiesAndForcesMT(tid);
    pthread_barrier_wait(&barrier);
    ComputeDensitiesMT(tid);
    pthread_barrier_wait(&barrier);
    ComputeDensities2MT(tid);
    pthread_barrier_wait(&barrier);
    ComputeForcesMT(tid);
    pthread_barrier_wait(&barrier);
    ProcessCollisionsMT(tid);
    pthread_barrier_wait(&barrier);
    AdvanceParticlesMT(tid);
    pthread_barrier_wait(&barrier);
#if defined(USE_ImpeneratableWall)
    // N.B. The integration of the position can place the particle
    // outside the domain. We now make a pass on the perimiter cells
    // to account for particle migration beyond domain.
    ProcessCollisions2MT(tid);
    pthread_barrier_wait(&barrier);
#endif
}

#ifndef ENABLE_VISUALIZATION
void *AdvanceFramesMT(void *args) {
    thread_args *targs = (thread_args *)args;
    std::cout << "tid: " << targs->tid << ", frames: " << targs->frames << std::endl;

    for (int i = 0; i < targs->frames; ++i) {
        AdvanceFrameMT(targs->tid);
    }

    return NULL;
}
#else
// Frame advancement function for worker threads
void *AdvanceFramesMT(void *args) {
    thread_args *targs = (thread_args *)args;

#if 1
    while (1)
#else
    for (int i = 0; i < targs->frames; ++i)
#endif
    {
        pthread_barrier_wait(&visualization_barrier);
        // Phase 1: Compute frame, visualization code blocked
        AdvanceFrameMT(targs->tid);
        pthread_barrier_wait(&visualization_barrier);
        // Phase 2: Visualize, worker threads blocked
    }

    return NULL;
}

// Frame advancement function for master thread (executes serial visualization code)
void AdvanceFrameVisualization() {
    // End of phase 2: Worker threads blocked, visualization code busy (last frame)
    pthread_barrier_wait(&visualization_barrier);
    // Phase 1: Visualization thread blocked, worker threads busy (next frame)
    pthread_barrier_wait(&visualization_barrier);
    // Begin of phase 2: Worker threads blocked, visualization code busy (next frame)
}
#endif  // ENABLE_VISUALIZATION

////////////////////////////////////////////////////////////////////////////////

std::vector<int> findChasOfCell(const Cell* cells)
{
    std::vector<int> res;

    const double* inspector = reinterpret_cast<const double*>(cells);
    for(int i = 0; i < sizeof(Cell) / sizeof(double); i += (CACHELINE_SIZE / sizeof(double))) {
        const auto cha = findCHAByHashing(reinterpret_cast<uintptr_t>(&inspector[i]));
        // std::cout << "i: " << i << ", cha: " << cha << std::endl;
        res.push_back(cha);
    }

    return res;
}

int main(int argc, char *argv[]) {
    using namespace std;
    // Cell is of size 896 bytes: 14 whole cache lines. this means that a single Cell might potentially have its coherence managed by 14 different CHAs.
    std::cout << "Cell size: " << sizeof(Cell) << std::endl;

    int res = posix_memalign((void **)(&cells), CACHELINE_SIZE, sizeof(struct Cell) * 1);
    if(res == 0) {
        std::cout << "posix_memalign success." << std::endl;
    } else {
        std::cerr << "posix_memalign errror.\n";
    }

    assert(static_cast<Cell*>(cells) != nullptr);


    const auto cell_chas = findChasOfCell(cells);
    // Cell is of size 896 bytes: 14 whole cache lines. this means that a single Cell will have its coherence managed by 14 CHAs. (not guaranteed to be distinct CHAs though)
    assert(cell_chas.size() == 14);
    for(const auto cell_cha : cell_chas) {
        std::cout << "cha: " << cell_cha << std::endl;
    }

    // assertRoot();  // no need for this to find chas by hashing.

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
    std::cout << "PARSEC Benchmark Suite Version "__PARSEC_XSTRING(PARSEC_VERSION) << std::endl << std::flush;
#else
    std::cout << "PARSEC Benchmark Suite" << std::endl << std::flush;
#endif  // PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
    __parsec_bench_begin(__parsec_fluidanimate);
#endif

    if (argc < 4 || argc >= 6) {
        std::cout << "Usage: " << argv[0] << " <threadnum> <framenum> <.fluid input file> [.fluid output file]"
                  << std::endl;
        return -1;
    }

    int threadnum = atoi(argv[1]);
    int framenum = atoi(argv[2]);
    std::cout << "threadnum: " << threadnum << ", framenum: " << framenum << std::endl;

    // Check arguments
    if (threadnum < 1) {
        std::cerr << "<threadnum> must at least be 1" << std::endl;
        return -1;
    }
    if (framenum < 1) {
        std::cerr << "<framenum> must at least be 1" << std::endl;
        return -1;
    }

#ifdef ENABLE_CFL_CHECK
    std::cout
        << "WARNING: Check for Courant–Friedrichs–Lewy condition enabled. Do not use for performance measurements."
        << std::endl;
#endif

    InitSim(argv[3], threadnum);
#ifdef ENABLE_VISUALIZATION
    InitVisualizationMode(&argc, argv, &AdvanceFrameVisualization, &numCells, &cells, &cnumPars);
#endif

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_begin();
#endif
#if defined(WIN32)
    thread_args *targs = (thread_args *)alloca(sizeof(thread_args) * threadnum);
#else
    thread_args targs[threadnum];
#endif
    for (int i = 0; i < threadnum; ++i) {
        targs[i].tid = i;
        targs[i].frames = framenum;
        pthread_create(&thread[i], &attr, AdvanceFramesMT, &targs[i]);
    }

    // *** PARALLEL PHASE *** //
#ifdef ENABLE_VISUALIZATION
    Visualize();
#endif

    for (int i = 0; i < threadnum; ++i) {
        pthread_join(thread[i], NULL);
    }
#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_end();
#endif

    if (argc > 4) SaveFile(argv[4]);
    CleanUpSim();

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_bench_end();
#endif
    // for (const auto &[threadid, addresses] : threadid_addresses_map) {
    //     std::cout << "thread id: " << threadid << " -- ";
    //     for (const auto &address : addresses) {
    //         std::cout << address << " -- ";
    //     }
    //     std::cout << std::endl;
    // }
    std::cout << "total thread count: " << threadid_addresses_map.size() << std::endl;

    assert(threadnum > 1);  // below algo depends on this.
    auto head = threadid_addresses_map.begin();
    auto tail = std::next(threadid_addresses_map.begin());

    // this is ranked_communication_count_per_pair wrt spmv repo.
    multiset<tuple<int, int, int>, greater<>>
        total_comm_count_t1_t2;  // set should suffice (compared to multiset). no tuple will be
                                 // the same since thread pairs are unique at this point here. but now, will make it multiset
    std::multiset<tuple<int, int, int, int>, greater<>> total_cha_freq_count_t1_t2;

    map<pair<int, int>, multiset<Cell *>> pairing_addresses;
    while (head != threadid_addresses_map.end()) {
        const auto orig_tail = tail;
        while (tail != threadid_addresses_map.end()) {
            const int t1 = head->first;
            const int t2 = tail->first;
            cout << "head: " << t1 << ", tail: " << t2 << endl;

            const multiset<Cell *> t1_addresses = head->second;
            const multiset<Cell *> t2_addresses = tail->second;

            std::multiset<Cell *> common_addresses;
            std::set_intersection(t1_addresses.begin(), t1_addresses.end(), t2_addresses.begin(), t2_addresses.end(),
                                  std::inserter(common_addresses, common_addresses.begin()));

            std::unordered_map<int, int> cha_freq_map;
            for(const Cell* common_addr : common_addresses) {
                for(const auto cha : findChasOfCell(common_addr)) {
                    ++cha_freq_map[cha];
                }
            }
            for(const auto& [cha, freq] : cha_freq_map) {
                total_cha_freq_count_t1_t2.insert({freq, cha, t1, t2});
            }

            total_comm_count_t1_t2.insert({common_addresses.size(), t1, t2});
            
            pairing_addresses[{t1, t2}] = common_addresses; // pairing is not used at the moment. here just for clarity.
            ++tail;
        }
        tail = std::next(orig_tail);
        ++head;
    }

    // for (const auto &[thread_pairs, common_addresses] : pairing_addresses) {
    //     const auto t1 = thread_pairs.first;
    //     const auto t2 = thread_pairs.second;
    //     const auto common_address_count = common_addresses.size();
    //     std::cout << "threads " << t1 << " and " << t2 << " have " << common_address_count << " common addresses"
    //               << endl;
    //     total_comm_count_t1_t2.insert({common_address_count, t1, t2});
    // }

    for (const auto &[total_comm_count, t1, t2] : total_comm_count_t1_t2) {
        std::cout << "total comm count: " << total_comm_count << ", t1: " << t1 << ", t2: " << t2 << endl;
    }
    for (const auto &[freq, cha, t1, t2] : total_cha_freq_count_t1_t2) {
        std::cout << "freq: " << freq << ", cha: " << cha << ", t1: " << t1 << ", t2: " << t2 << endl;
    }    

    int mapped_thread_count = 0;
    auto it = total_cha_freq_count_t1_t2.begin();
    auto it1 = total_comm_count_t1_t2.begin();
    std::vector<int> thread_to_core(threadnum, -1);

    // fprintf(stderr, "before topology creation\n");
    auto topo = Topology(cha_core_map, CAPID6);
    topo.printTopology();
    std::vector<Tile> mapped_tiles;
    // SPDLOG_TRACE("~~~~~~~~~~~~~~~~");
    //  fprintf(stderr, "before thread mapping creation\n");

    // start
    // it = ranked_cha_access_count_per_pair.begin();
    while (mapped_tiles.size() < threadnum &&
           /*it != ranked_cha_access_count_per_pair.end()*/ it1 != total_comm_count_t1_t2.end()) {
        // std::pair<int, int> tid_pair(std::get<2>(*it), std::get<3>(*it));
        std::pair<int, int> tid_pair(std::get<1>(*it1), std::get<2>(*it1));
        if (thread_to_core[tid_pair.first] == -1 && thread_to_core[tid_pair.second] == -1) {
            // SPDLOG_TRACE("cha with max access: {}", std::get<1>(*it));
            int cha_id = getMostAccessedCHA(tid_pair.first, tid_pair.second, total_cha_freq_count_t1_t2, topo);
            if (cha_id == -1) {
                // SPDLOG_INFO("error: cha is -1");
                it1++;
                continue;
            }
            // auto tile = topo.getTile(std::get<1>(*it));
            auto tile = topo.getTile(cha_id);
            // SPDLOG_TRACE("cha {}, is colocated with core {}", cha_id, tile.core);
            // if (thread_to_core[tid_pair.first] == -1)
            {
                // SPDLOG_INFO("fetching a tile closest to tile with cha {} and core {}, cha supposed to be {}",
                // tile.cha, tile.core, std::get<1>(*it));
                auto closest_tile = topo.getClosestTile(tile, mapped_tiles);
                // auto closest_tile = topo.getClosestTilewithThreshold(tile, mapped_tiles);
                // SPDLOG_TRACE("* closest _available_ core to cha {} is: {}", tile.cha, closest_tile.core);
                mapped_tiles.push_back(closest_tile);
                thread_to_core[tid_pair.first] = closest_tile.core;
                // SPDLOG_TRACE("assigned thread with id {} to core {}", tid_pair.first, closest_tile.core);
            }
#if 0
            else
            {
                SPDLOG_TRACE("--> Already assigned thread with id {} to core {}, skipping it.", tid_pair.first, thread_to_core[tid_pair.first]);
            }
#endif

            // if (thread_to_core[tid_pair.second] == -1)
            {
                // SPDLOG_INFO("fetching a tile closest to tile with cha {} and core {}, cha supposed to be {}",
                // tile.cha, tile.core, std::get<1>(*it));
                auto closest_tile = topo.getClosestTile(tile, mapped_tiles);
                // auto closest_tile = topo.getClosestTilewithThreshold(tile, mapped_tiles);
                // SPDLOG_TRACE("# closest _available_ core to cha {} is: {}", tile.cha, closest_tile.core);
                mapped_tiles.push_back(closest_tile);
                thread_to_core[tid_pair.second] = closest_tile.core;
                // SPDLOG_TRACE("assigned thread with id {} to core {}", tid_pair.second, closest_tile.core);
            }
#if 0
            else
            {
                SPDLOG_TRACE("--> Already assigned thread with id {} to core {}, skipping it.", tid_pair.second, thread_to_core[tid_pair.second]);
            }
#endif
        }
        //#if 0
        else if (thread_to_core[tid_pair.first] == -1) {
            auto tile = topo.getTileByCore(thread_to_core[tid_pair.second]);
            auto closest_tile = topo.getClosestTile(tile, mapped_tiles);
            mapped_tiles.push_back(closest_tile);
            thread_to_core[tid_pair.first] = closest_tile.core;
        } else if (thread_to_core[tid_pair.second] == -1) {
            auto tile = topo.getTileByCore(thread_to_core[tid_pair.first]);
            auto closest_tile = topo.getClosestTile(tile, mapped_tiles);
            mapped_tiles.push_back(closest_tile);
            thread_to_core[tid_pair.second] = closest_tile.core;
        }
        //#endif

        it1++;
    }
    // end
    int i = 0;
    for (auto ptr : thread_to_core) {
        std::cout << "thread " << i << " is mapped to core " << ptr << std::endl;
        // SPDLOG_INFO("thread {} is mapped to core {} ", i, ptr);
        i++;
    }

    // -> thread to core is the optimized binding. thread-0 will be bound to core thread_to_core[0] etc.

    //////////////////// TODOs /////////////////

    // TODO: using "pairing_addresses", create ranked_cha_access_count_per_pair
    // set<tuple<cha_freq, cha, tid-a, tid-b>>

    // TODO: copy the mapping algo here.

    // TODO: warm cache before the benchmark so both versions can benefit from cache access initially.

    // TODO: measure preprocessing overhead.

    // TODO: compile on koc cascade.

    // TODO: run in_500K.fluid and in_300K.fluid files and benchmark it!

    // TODO: compare "optimized" version with the default binding (just bind threads to cores 0 to n - 1, n being the
    // total core count).

    // TODO: write a bash script and let the app run multiple times (at least 100). each run should last longer than 5
    // seconds.

    // TODO: I should link to numa library and allocate memory via numa_alloc_onnode???
    // I might not need this though: In Linux, the default behavior for memory allocation is "first-touch". This means
    // that the memory will not be actually allocated to a specific NUMA node until the application first accesses it
    // (e.g., writes to it). At that point, it will be allocated on the NUMA node local to the core that did the
    // accessing. So, if core 0 (which is on NUMA node 0) is the first to access the memory, then it will be allocated
    // on NUMA node 0.

    return 0;

    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();
    // long_operation();
    auto t2 = high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = duration_cast<milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_int.count() << "ms\n";
    std::cout << ms_double.count() << "ms\n";
}

////////////////////////////////////////////////////////////////////////////////
