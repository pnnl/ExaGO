#include <algorithm>
#include <cmath>
#include <iostream>

#include <mpi.h>

#include <psimpl.h>
#include <psse.hpp>

int TEST_RESULT = EXIT_SUCCESS;

#define TEST(cond)                                                             \
  {                                                                            \
    if (!(cond)) {                                                             \
      __local.result = false;                                                  \
      TEST_RESULT = EXIT_FAILURE;                                              \
      std::cout << "Test failed: line " << __LINE__ << '\n';                   \
    }                                                                          \
  }

#define DEFAULT_RTOL 1.0e-5

template <typename T> bool CheckClose(T a, T b, T rtol = DEFAULT_RTOL) {
  return std::fabs(a - b) <= rtol * std::max(std::fabs(a), std::fabs(b));
}

#define TEST_CLOSE_FAIL_MSG(a, b, ...)                                         \
  std::cout << "\tGot: " << a << " == " << b << '\n';

#define TEST_CLOSE(...)                                                        \
  [&]() {                                                                      \
    bool __chk = CheckClose(__VA_ARGS__);                                      \
    TEST(__chk);                                                               \
    if (!__chk) {                                                              \
      TEST_CLOSE_FAIL_MSG(__VA_ARGS__);                                        \
    }                                                                          \
  }()

#define TEST_EQUAL(a, b)                                                       \
  [&]() {                                                                      \
    bool __chk = (a == b);                                                     \
    TEST(__chk);                                                               \
    if (!__chk) {                                                              \
      TEST_CLOSE_FAIL_MSG(a, b);                                               \
    }                                                                          \
  }()

std::string Strip(std::string str) {
  auto notspace = [](char c) { return !std::isspace(c); };
  str.erase(begin(str), std::find_if(begin(str), end(str), notspace));
  str.erase(std::find_if(rbegin(str), rend(str), notspace).base(), end(str));
  return str;
}

auto Deg2Rad = [](auto &&deg) { return deg * M_PI / 180.0; };

struct PSHolder {
  PS ps;

  PSHolder() { PSCreate(MPI_COMM_WORLD, &ps); }

  ~PSHolder() { PSDestroy(&ps); }
};

PSHolder ReadPSData(const std::string &filename) {
  PSHolder psh;
  PSReadPSSERawData(psh.ps, filename.c_str());
  return psh;
}

PSHolder NetworkToPS(const exago::psse::Network &nw) {
  PSHolder psh;
  exago::psse::ConvertToPS(psh.ps, nw);
  return psh;
}

#define CHECK(cond)                                                            \
  {                                                                            \
    if (!(cond)) {                                                             \
      __local.result = false;                                                  \
      std::cout << "Test failed: line " << __LINE__ << '\n';                   \
    }                                                                          \
  }

struct LocalResult {
  mutable bool result{true};
};

#define TEST_FUNCTION(name) auto name = [__local = LocalResult()]

#define TEST_FUNCTION_RETURN                                                   \
  bool __result = __local.result;                                              \
  __local.result = true;                                                       \
  return __result;
#define END_TEST_FUNCTION ;

TEST_FUNCTION(check_ps_data_ieee9bus)(PS ps) {
  auto MVAbase = ps->MVAbase;
  TEST_EQUAL(ps->MVAbase, 100.0);
  TEST_EQUAL(ps->nbus, 9);
  TEST_EQUAL(ps->maxbusnum, 9);

  // buses
  TEST_EQUAL(ps->Nbus, 9);

  TEST_EQUAL(ps->bus[0].bus_i, 1);
  TEST_EQUAL(Strip(ps->bus[0].name), "BUS1");
  TEST_EQUAL(ps->bus[0].basekV, 16.5);
  TEST_EQUAL(ps->bus[0].ide, 3);
  TEST_EQUAL(ps->bus[0].area, 1);
  TEST_EQUAL(ps->bus[0].zone, 1);
  TEST_EQUAL(ps->bus[0].owner, 1);
  TEST_EQUAL(ps->bus[0].vm, 1.04);
  TEST_CLOSE(ps->bus[0].va, 0.0);

  TEST_EQUAL(ps->bus[1].bus_i, 2);
  TEST_EQUAL(Strip(ps->bus[1].name), "BUS2");
  TEST_EQUAL(ps->bus[1].basekV, 18.0);
  TEST_EQUAL(ps->bus[1].ide, 2);
  TEST_EQUAL(ps->bus[1].area, 2);
  TEST_EQUAL(ps->bus[1].zone, 2);
  TEST_EQUAL(ps->bus[1].owner, 1);
  TEST_EQUAL(ps->bus[1].vm, 1.025);
  TEST_CLOSE(ps->bus[1].va, Deg2Rad(9.28));

  TEST_EQUAL(ps->bus[2].bus_i, 3);
  TEST_EQUAL(Strip(ps->bus[2].name), "BUS3");
  TEST_EQUAL(ps->bus[2].basekV, 13.8);
  TEST_EQUAL(ps->bus[2].ide, 2);
  TEST_EQUAL(ps->bus[2].area, 2);
  TEST_EQUAL(ps->bus[2].zone, 3);
  TEST_EQUAL(ps->bus[2].owner, 1);
  TEST_EQUAL(ps->bus[2].vm, 1.025);
  TEST_CLOSE(ps->bus[2].va, Deg2Rad(4.6648));

  TEST_EQUAL(ps->bus[3].bus_i, 4);
  TEST_EQUAL(Strip(ps->bus[3].name), "BUS4");
  TEST_EQUAL(ps->bus[3].basekV, 230.0);
  TEST_EQUAL(ps->bus[3].ide, 1);
  TEST_EQUAL(ps->bus[3].area, 1);
  TEST_EQUAL(ps->bus[3].zone, 1);
  TEST_EQUAL(ps->bus[3].owner, 1);
  TEST_EQUAL(ps->bus[3].vm, 1.02579);
  TEST_CLOSE(ps->bus[3].va, Deg2Rad(-2.2168));

  TEST_EQUAL(ps->bus[4].bus_i, 5);
  TEST_EQUAL(Strip(ps->bus[4].name), "BUS5");
  TEST_EQUAL(ps->bus[4].basekV, 230.0);
  TEST_EQUAL(ps->bus[4].ide, 1);
  TEST_EQUAL(ps->bus[4].area, 1);
  TEST_EQUAL(ps->bus[4].zone, 4);
  TEST_EQUAL(ps->bus[4].owner, 2);
  TEST_EQUAL(ps->bus[4].vm, 0.99563);
  TEST_CLOSE(ps->bus[4].va, Deg2Rad(-3.9888));

  TEST_EQUAL(ps->bus[5].bus_i, 6);
  TEST_EQUAL(Strip(ps->bus[5].name), "BUS6");
  TEST_EQUAL(ps->bus[5].basekV, 230.0);
  TEST_EQUAL(ps->bus[5].ide, 1);
  TEST_EQUAL(ps->bus[5].area, 1);
  TEST_EQUAL(ps->bus[5].zone, 5);
  TEST_EQUAL(ps->bus[5].owner, 2);
  TEST_EQUAL(ps->bus[5].vm, 1.01265);
  TEST_CLOSE(ps->bus[5].va, Deg2Rad(-3.6874));

  TEST_EQUAL(ps->bus[6].bus_i, 7);
  TEST_EQUAL(Strip(ps->bus[6].name), "BUS7");
  TEST_EQUAL(ps->bus[6].basekV, 230.0);
  TEST_EQUAL(ps->bus[6].ide, 1);
  TEST_EQUAL(ps->bus[6].area, 2);
  TEST_EQUAL(ps->bus[6].zone, 2);
  TEST_EQUAL(ps->bus[6].owner, 1);
  TEST_EQUAL(ps->bus[6].vm, 1.02577);
  TEST_CLOSE(ps->bus[6].va, Deg2Rad(3.7197));

  TEST_EQUAL(ps->bus[7].bus_i, 8);
  TEST_EQUAL(Strip(ps->bus[7].name), "BUS8");
  TEST_EQUAL(ps->bus[7].basekV, 230.0);
  TEST_EQUAL(ps->bus[7].ide, 1);
  TEST_EQUAL(ps->bus[7].area, 2);
  TEST_EQUAL(ps->bus[7].zone, 6);
  TEST_EQUAL(ps->bus[7].owner, 2);
  TEST_EQUAL(ps->bus[7].vm, 1.01588);
  TEST_CLOSE(ps->bus[7].va, Deg2Rad(0.7275));

  TEST_EQUAL(ps->bus[8].bus_i, 9);
  TEST_EQUAL(Strip(ps->bus[8].name), "BUS9");
  TEST_EQUAL(ps->bus[8].basekV, 230.0);
  TEST_EQUAL(ps->bus[8].ide, 1);
  TEST_EQUAL(ps->bus[8].area, 2);
  TEST_EQUAL(ps->bus[8].zone, 3);
  TEST_EQUAL(ps->bus[8].owner, 1);
  TEST_EQUAL(ps->bus[8].vm, 1.03235);
  TEST_CLOSE(ps->bus[8].va, Deg2Rad(1.9667));

  // loads
  TEST_EQUAL(ps->Nload, 3);

  TEST_EQUAL(ps->load[0].bus_i, 5);
  TEST_EQUAL(Strip(ps->load[0].id), "1");
  TEST_EQUAL(ps->load[0].status, 1);
  TEST_EQUAL(ps->load[0].area, 1);
  TEST_EQUAL(ps->load[0].zone, 1);
  TEST_CLOSE(ps->load[0].pl, 125.0 / MVAbase);
  TEST_CLOSE(ps->load[0].ql, 50.0 / MVAbase);
  TEST_EQUAL(ps->load[0].ip, 0.0);
  TEST_EQUAL(ps->load[0].iq, 0.0);
  TEST_EQUAL(ps->load[0].yp, 0.0);
  TEST_EQUAL(ps->load[0].yq, 0.0);
  TEST_EQUAL(ps->load[0].owner, 1);
  TEST_EQUAL(ps->load[0].scale, 1);

  TEST_EQUAL(ps->load[1].bus_i, 6);
  TEST_EQUAL(Strip(ps->load[1].id), "1");
  TEST_EQUAL(ps->load[1].status, 1);
  TEST_EQUAL(ps->load[1].area, 1);
  TEST_EQUAL(ps->load[1].zone, 1);
  TEST_CLOSE(ps->load[1].pl, 90.0 / MVAbase);
  TEST_CLOSE(ps->load[1].ql, 30.0 / MVAbase);
  TEST_EQUAL(ps->load[1].ip, 0.0);
  TEST_EQUAL(ps->load[1].iq, 0.0);
  TEST_EQUAL(ps->load[1].yp, 0.0);
  TEST_EQUAL(ps->load[1].yq, 0.0);
  TEST_EQUAL(ps->load[1].owner, 1);
  TEST_EQUAL(ps->load[1].scale, 1);

  TEST_EQUAL(ps->load[2].bus_i, 8);
  TEST_EQUAL(Strip(ps->load[2].id), "1");
  TEST_EQUAL(ps->load[2].status, 1);
  TEST_EQUAL(ps->load[2].area, 1);
  TEST_EQUAL(ps->load[2].zone, 1);
  TEST_CLOSE(ps->load[2].pl, 100.0 / MVAbase);
  TEST_CLOSE(ps->load[2].ql, 35.0 / MVAbase);
  TEST_EQUAL(ps->load[2].ip, 0.0);
  TEST_EQUAL(ps->load[2].iq, 0.0);
  TEST_EQUAL(ps->load[2].yp, 0.0);
  TEST_EQUAL(ps->load[2].yq, 0.0);
  TEST_EQUAL(ps->load[2].owner, 1);
  TEST_EQUAL(ps->load[2].scale, 1);

  // generators
  TEST_EQUAL(ps->Ngen, 3);

  TEST_EQUAL(ps->gen[0].bus_i, 1);
  TEST_EQUAL(Strip(ps->gen[0].id), "1");
  TEST_CLOSE(ps->gen[0].pg, 71.641 / MVAbase);
  TEST_CLOSE(ps->gen[0].qg, 27.046 / MVAbase);
  TEST_CLOSE(ps->gen[0].qt, 300.0 / MVAbase);
  TEST_CLOSE(ps->gen[0].qb, -300.0 / MVAbase);
  TEST_EQUAL(ps->gen[0].vs, 1.04);
  TEST_EQUAL(ps->gen[0].ireg, 0);
  TEST_EQUAL(ps->gen[0].mbase, 260.0);
  TEST_EQUAL(ps->gen[0].zr, 1.0e-4);
  TEST_EQUAL(ps->gen[0].zx, 0.1);
  TEST_EQUAL(ps->gen[0].rt, 0.0);
  TEST_EQUAL(ps->gen[0].xt, 0.0);
  TEST_EQUAL(ps->gen[0].gtap, 1.0);
  TEST_EQUAL(ps->gen[0].status, 1);
  TEST_EQUAL(ps->gen[0].rmpct, 100.0);
  TEST_CLOSE(ps->gen[0].pt, 250.0 / MVAbase);
  TEST_CLOSE(ps->gen[0].pb, 10.0 / MVAbase);
  TEST_EQUAL(ps->gen[0].o1, 1);
  TEST_EQUAL(ps->gen[0].f1, 1.0);

  TEST_EQUAL(ps->gen[1].bus_i, 2);
  TEST_EQUAL(Strip(ps->gen[1].id), "1");
  TEST_CLOSE(ps->gen[1].pg, 163.0 / MVAbase);
  TEST_CLOSE(ps->gen[1].qg, 6.654 / MVAbase);
  TEST_CLOSE(ps->gen[1].qt, 300.0 / MVAbase);
  TEST_CLOSE(ps->gen[1].qb, -300.0 / MVAbase);
  TEST_EQUAL(ps->gen[1].vs, 1.025);
  TEST_EQUAL(ps->gen[1].ireg, 0);
  TEST_EQUAL(ps->gen[1].mbase, 310.0);
  TEST_EQUAL(ps->gen[1].zr, 1.0e-4);
  TEST_EQUAL(ps->gen[1].zx, 0.21);
  TEST_EQUAL(ps->gen[1].rt, 0.0);
  TEST_EQUAL(ps->gen[1].xt, 0.0);
  TEST_EQUAL(ps->gen[1].gtap, 1.0);
  TEST_EQUAL(ps->gen[1].status, 1);
  TEST_EQUAL(ps->gen[1].rmpct, 100.0);
  TEST_CLOSE(ps->gen[1].pt, 300.0 / MVAbase);
  TEST_CLOSE(ps->gen[1].pb, 10.0 / MVAbase);
  TEST_EQUAL(ps->gen[1].o1, 1);
  TEST_EQUAL(ps->gen[1].f1, 1.0);

  TEST_EQUAL(ps->gen[2].bus_i, 3);
  TEST_EQUAL(Strip(ps->gen[2].id), "1");
  TEST_CLOSE(ps->gen[2].pg, 85.0 / MVAbase);
  TEST_CLOSE(ps->gen[2].qg, -10.86 / MVAbase);
  TEST_CLOSE(ps->gen[2].qt, 300.0 / MVAbase);
  TEST_CLOSE(ps->gen[2].qb, -300.0 / MVAbase);
  TEST_EQUAL(ps->gen[2].vs, 1.025);
  TEST_EQUAL(ps->gen[2].ireg, 0);
  TEST_EQUAL(ps->gen[2].mbase, 280.0);
  TEST_EQUAL(ps->gen[2].zr, 1.0e-4);
  TEST_EQUAL(ps->gen[2].zx, 0.21);
  TEST_EQUAL(ps->gen[2].rt, 0.0);
  TEST_EQUAL(ps->gen[2].xt, 0.0);
  TEST_EQUAL(ps->gen[2].gtap, 1.0);
  TEST_EQUAL(ps->gen[2].status, 1);
  TEST_EQUAL(ps->gen[2].rmpct, 100.0);
  TEST_CLOSE(ps->gen[2].pt, 270.0 / MVAbase);
  TEST_CLOSE(ps->gen[2].pb, 10.0 / MVAbase);
  TEST_EQUAL(ps->gen[2].o1, 1);
  TEST_EQUAL(ps->gen[2].f1, 1.0);

  // lines
  TEST_EQUAL(ps->Nline, 9);

  // branches
  TEST_EQUAL(ps->line[0].fbus, 4);
  TEST_EQUAL(ps->line[0].tbus, 5);
  TEST_EQUAL(Strip(ps->line[0].ckt), "1");
  TEST_EQUAL(ps->line[0].r, 0.01);
  TEST_EQUAL(ps->line[0].x, 0.085);
  TEST_EQUAL(ps->line[0].b, 0.176);
  TEST_EQUAL(ps->line[0].rateA, 250.0);
  TEST_EQUAL(ps->line[0].rateB, 250.0);
  TEST_EQUAL(ps->line[0].rateC, 250.0);
  TEST_EQUAL(ps->line[0].gi, 0.0);
  TEST_EQUAL(ps->line[0].bi, 0.0);
  TEST_EQUAL(ps->line[0].gj, 0.0);
  TEST_EQUAL(ps->line[0].bj, 0.0);
  TEST_EQUAL(ps->line[0].status, 1);
  // MET not retained
  TEST_EQUAL(ps->line[0].length, 0.0);
  TEST_EQUAL(ps->line[0].o1, 1);
  TEST_EQUAL(ps->line[0].f1, 1.0);

  TEST_EQUAL(ps->line[1].fbus, 4);
  TEST_EQUAL(ps->line[1].tbus, 6);
  TEST_EQUAL(Strip(ps->line[1].ckt), "1");
  TEST_EQUAL(ps->line[1].r, 0.017);
  TEST_EQUAL(ps->line[1].x, 0.092);
  TEST_EQUAL(ps->line[1].b, 0.158);
  TEST_EQUAL(ps->line[1].rateA, 250.0);
  TEST_EQUAL(ps->line[1].rateB, 250.0);
  TEST_EQUAL(ps->line[1].rateC, 250.0);
  TEST_EQUAL(ps->line[1].gi, 0.0);
  TEST_EQUAL(ps->line[1].bi, 0.0);
  TEST_EQUAL(ps->line[1].gj, 0.0);
  TEST_EQUAL(ps->line[1].bj, 0.0);
  TEST_EQUAL(ps->line[1].status, 1);
  // MET not retained
  TEST_EQUAL(ps->line[1].length, 0.0);
  TEST_EQUAL(ps->line[1].o1, 1);
  TEST_EQUAL(ps->line[1].f1, 1.0);

  TEST_EQUAL(ps->line[2].fbus, 5);
  TEST_EQUAL(ps->line[2].tbus, 7);
  TEST_EQUAL(Strip(ps->line[2].ckt), "1");
  TEST_EQUAL(ps->line[2].r, 0.032);
  TEST_EQUAL(ps->line[2].x, 0.161);
  TEST_EQUAL(ps->line[2].b, 0.306);
  TEST_EQUAL(ps->line[2].rateA, 250.0);
  TEST_EQUAL(ps->line[2].rateB, 250.0);
  TEST_EQUAL(ps->line[2].rateC, 250.0);
  TEST_EQUAL(ps->line[2].gi, 0.0);
  TEST_EQUAL(ps->line[2].bi, 0.0);
  TEST_EQUAL(ps->line[2].gj, 0.0);
  TEST_EQUAL(ps->line[2].bj, 0.0);
  TEST_EQUAL(ps->line[2].status, 1);
  // MET not retained
  TEST_EQUAL(ps->line[2].length, 0.0);
  TEST_EQUAL(ps->line[2].o1, 1);
  TEST_EQUAL(ps->line[2].f1, 1.0);

  TEST_EQUAL(ps->line[3].fbus, 6);
  TEST_EQUAL(ps->line[3].tbus, 9);
  TEST_EQUAL(Strip(ps->line[3].ckt), "1");
  TEST_EQUAL(ps->line[3].r, 0.039);
  TEST_EQUAL(ps->line[3].x, 0.17);
  TEST_EQUAL(ps->line[3].b, 0.358);
  TEST_EQUAL(ps->line[3].rateA, 150.0);
  TEST_EQUAL(ps->line[3].rateB, 150.0);
  TEST_EQUAL(ps->line[3].rateC, 150.0);
  TEST_EQUAL(ps->line[3].gi, 0.0);
  TEST_EQUAL(ps->line[3].bi, 0.0);
  TEST_EQUAL(ps->line[3].gj, 0.0);
  TEST_EQUAL(ps->line[3].bj, 0.0);
  TEST_EQUAL(ps->line[3].status, 1);
  // MET not retained
  TEST_EQUAL(ps->line[3].length, 0.0);
  TEST_EQUAL(ps->line[3].o1, 1);
  TEST_EQUAL(ps->line[3].f1, 1.0);

  TEST_EQUAL(ps->line[4].fbus, 7);
  TEST_EQUAL(ps->line[4].tbus, 8);
  TEST_EQUAL(Strip(ps->line[4].ckt), "1");
  TEST_EQUAL(ps->line[4].r, 0.0085);
  TEST_EQUAL(ps->line[4].x, 0.072);
  TEST_EQUAL(ps->line[4].b, 0.149);
  TEST_EQUAL(ps->line[4].rateA, 250.0);
  TEST_EQUAL(ps->line[4].rateB, 250.0);
  TEST_EQUAL(ps->line[4].rateC, 250.0);
  TEST_EQUAL(ps->line[4].gi, 0.0);
  TEST_EQUAL(ps->line[4].bi, 0.0);
  TEST_EQUAL(ps->line[4].gj, 0.0);
  TEST_EQUAL(ps->line[4].bj, 0.0);
  TEST_EQUAL(ps->line[4].status, 1);
  // MET not retained
  TEST_EQUAL(ps->line[4].length, 0.0);
  TEST_EQUAL(ps->line[4].o1, 1);
  TEST_EQUAL(ps->line[4].f1, 1.0);

  TEST_EQUAL(ps->line[5].fbus, 8);
  TEST_EQUAL(ps->line[5].tbus, 9);
  TEST_EQUAL(Strip(ps->line[5].ckt), "1");
  TEST_EQUAL(ps->line[5].r, 0.0119);
  TEST_EQUAL(ps->line[5].x, 0.1008);
  TEST_EQUAL(ps->line[5].b, 0.209);
  TEST_EQUAL(ps->line[5].rateA, 150.0);
  TEST_EQUAL(ps->line[5].rateB, 150.0);
  TEST_EQUAL(ps->line[5].rateC, 150.0);
  TEST_EQUAL(ps->line[5].gi, 0.0);
  TEST_EQUAL(ps->line[5].bi, 0.0);
  TEST_EQUAL(ps->line[5].gj, 0.0);
  TEST_EQUAL(ps->line[5].bj, 0.0);
  TEST_EQUAL(ps->line[5].status, 1);
  // MET not retained
  TEST_EQUAL(ps->line[5].length, 0.0);
  TEST_EQUAL(ps->line[5].o1, 1);
  TEST_EQUAL(ps->line[5].f1, 1.0);

  // transformers
  TEST_EQUAL(ps->line[6].fbus, 1);
  TEST_EQUAL(ps->line[6].tbus, 4);
  // K not retained
  TEST_EQUAL(Strip(ps->line[6].ckt), "T1");
  // CW, CZ, CM, MAG1, MAG2, NMETR, NAME not retained
  TEST_EQUAL(ps->line[6].status, 1);
  TEST_EQUAL(ps->line[6].o1, 1);
  TEST_EQUAL(ps->line[6].f1, 1.0);
  // file line 2
  TEST_EQUAL(ps->line[6].r, 0.0);
  TEST_EQUAL(ps->line[6].x, 0.0576);
  TEST_EQUAL(ps->line[6].sbase12, 100.0);
  // file line 3
  TEST_EQUAL(ps->line[6].tapratio, 1.0);
  // NOMV1 not retained
  TEST_EQUAL(ps->line[6].phaseshift, 0.0);
  TEST_EQUAL(ps->line[6].rateA, 0.0);
  TEST_EQUAL(ps->line[6].rateB, 0.0);
  TEST_EQUAL(ps->line[6].rateC, 0.0);
  // COD1, CONT1, RMA1, RMI1, VMA1, VMI1, NTP1, TAB1, CR1, CX1, CNXA1 not
  // retained record 4 WINDV2, NOMV2 not retained

  TEST_EQUAL(ps->line[7].fbus, 2);
  TEST_EQUAL(ps->line[7].tbus, 7);
  // K not retained
  TEST_EQUAL(Strip(ps->line[7].ckt), "T2");
  // CW, CZ, CM, MAG1, MAG2, NMETR, NAME not retained
  TEST_EQUAL(ps->line[7].status, 1);
  TEST_EQUAL(ps->line[7].o1, 1);
  TEST_EQUAL(ps->line[7].f1, 1.0);
  // record 2
  TEST_EQUAL(ps->line[7].r, 0.0);
  TEST_EQUAL(ps->line[7].x, 0.0625);
  TEST_EQUAL(ps->line[7].sbase12, 100.0);
  // record 3
  TEST_EQUAL(ps->line[7].tapratio, 1.0);
  // NOMV1 not retained
  TEST_EQUAL(ps->line[7].phaseshift, 0.0);
  TEST_EQUAL(ps->line[7].rateA, 0.0);
  TEST_EQUAL(ps->line[7].rateB, 0.0);
  TEST_EQUAL(ps->line[7].rateC, 0.0);
  // COD1, CONT1, RMA1, RMI1, VMA1, VMI1, NTP1, TAB1, CR1, CX1, CNXA1 not
  // retained record 4 WINDV2, NOMV2 not retained

  TEST_EQUAL(ps->line[8].fbus, 9);
  TEST_EQUAL(ps->line[8].tbus, 3);
  // K not retained
  TEST_EQUAL(Strip(ps->line[8].ckt), "T3");
  // CW, CZ, CM, MAG1, MAG2, NMETR, NAME not retained
  TEST_EQUAL(ps->line[8].status, 1);
  TEST_EQUAL(ps->line[8].o1, 1);
  TEST_EQUAL(ps->line[8].f1, 1.0);
  // record 2
  TEST_EQUAL(ps->line[8].r, 0.0);
  TEST_EQUAL(ps->line[8].x, 0.0586);
  TEST_EQUAL(ps->line[8].sbase12, 100.0);
  // record 3
  TEST_EQUAL(ps->line[8].tapratio, 1.0);
  // NOMV1 not retained
  TEST_EQUAL(ps->line[8].phaseshift, 0.0);
  TEST_EQUAL(ps->line[8].rateA, 0.0);
  TEST_EQUAL(ps->line[8].rateB, 0.0);
  TEST_EQUAL(ps->line[8].rateC, 0.0);
  // COD1, CONT1, RMA1, RMI1, VMA1, VMI1, NTP1, TAB1, CR1, CX1, CNXA1 not
  // retained record 4 WINDV2, NOMV2 not retained

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(check_equal_bus)(PSBUS a, PSBUS b) {
  TEST_EQUAL(a->bus_i, b->bus_i);
  TEST_EQUAL(Strip(a->i), Strip(b->i));
  TEST_EQUAL(Strip(a->name), Strip(b->name));
  TEST_EQUAL(a->basekV, b->basekV);
  TEST_EQUAL(a->ide, b->ide);
  TEST_EQUAL(a->gl, b->gl);
  TEST_EQUAL(a->bl, b->bl);
  TEST_EQUAL(a->area, b->area);
  TEST_EQUAL(a->zone, b->zone);
  TEST_EQUAL(a->vm, b->vm);
  TEST_CLOSE(a->va, b->va);
  TEST_EQUAL(a->owner, b->owner);
  TEST_EQUAL(a->Vmax, b->Vmax);
  TEST_EQUAL(a->Vmin, b->Vmin);
  TEST_EQUAL(a->nvhi, b->nvhi);
  TEST_EQUAL(a->nvlo, b->nvlo);
  TEST_EQUAL(a->evhi, b->evhi);
  TEST_EQUAL(a->evlo, b->evlo);
  TEST_EQUAL(a->internal_i, b->internal_i);
  TEST_EQUAL(a->ngen, b->ngen);
  TEST_EQUAL(a->ngenON, b->ngenON);
  for (int i = 0; i < a->ngen; ++i) {
    TEST_EQUAL(a->gidx[i], b->gidx[i]);
  }
  TEST_EQUAL(a->nload, b->nload);
  for (int i = 0; i < a->nload; ++i) {
    TEST_EQUAL(a->lidx[i], b->lidx[i]);
  }
  TEST_EQUAL(a->nshunt, b->nshunt);
  TEST_EQUAL(a->qrange, b->qrange);
  TEST_EQUAL(a->qmintot, b->qmintot);
  TEST_EQUAL(a->Pgtot, b->Pgtot);
  TEST_EQUAL(a->MVAbasetot, b->MVAbasetot);
  TEST_EQUAL(a->nconnlines, b->nconnlines);
  for (int i = 0; i < a->nconnlines; ++i) {
    TEST_EQUAL(a->connlines[i], b->connlines[i]);
  }
  TEST_EQUAL(a->isghost, b->isghost);

  // gens, loads

  TEST_EQUAL(a->pimb, b->pimb);
  TEST_EQUAL(a->qimb, b->qimb);
  TEST_EQUAL(a->mult_pmis, b->mult_pmis);
  TEST_EQUAL(a->mult_qmis, b->mult_qmis);
  TEST_EQUAL(a->startloc, b->startloc);
  TEST_EQUAL(a->startlocglob, b->startlocglob);
  TEST_EQUAL(a->nxV, b->nxV);
  TEST_EQUAL(a->nxshunt, b->nxshunt);
  TEST_EQUAL(a->nxpimb, b->nxpimb);
  TEST_EQUAL(a->startxVloc, b->startxVloc);
  TEST_EQUAL(a->startxshuntloc, b->startxshuntloc);
  TEST_EQUAL(a->startxpimbloc, b->startxpimbloc);
  TEST_EQUAL(a->startxVlocglob, b->startxVlocglob);
  TEST_EQUAL(a->startxshuntlocglob, b->startxshuntlocglob);
  TEST_EQUAL(a->startxpimblocglob, b->startxpimblocglob);
  TEST_EQUAL(a->nx, b->nx);
  TEST_EQUAL(a->nconeq, b->nconeq);
  TEST_EQUAL(a->nconineq, b->nconineq);
  TEST_EQUAL(a->starteqloc, b->starteqloc);
  TEST_EQUAL(a->startineqloc, b->startineqloc);

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(check_equal_load)(PSLOAD a, PSLOAD b) {
  TEST_EQUAL(a->bus_i, b->bus_i);
  TEST_EQUAL(Strip(a->i), Strip(b->i));
  TEST_EQUAL(Strip(a->id), Strip(b->id));
  TEST_EQUAL(a->status, b->status);
  TEST_EQUAL(a->area, b->area);
  TEST_EQUAL(a->zone, b->zone);
  TEST_EQUAL(a->pl, b->pl);
  TEST_EQUAL(a->ql, b->ql);
  TEST_EQUAL(a->ip, b->ip);
  TEST_EQUAL(a->iq, b->iq);
  TEST_EQUAL(a->yp, b->yp);
  TEST_EQUAL(a->yq, b->yq);
  TEST_EQUAL(a->owner, b->owner);
  TEST_EQUAL(a->internal_i, b->internal_i);
  TEST_EQUAL(a->scale, b->scale);
  TEST_EQUAL(a->intrpt, b->intrpt);
  TEST_EQUAL(a->pl_loss, b->pl_loss);
  TEST_EQUAL(a->ql_loss, b->ql_loss);
  TEST_EQUAL(a->loss_cost, b->loss_cost);
  TEST_EQUAL(a->loss_frac, b->loss_frac);
  TEST_EQUAL(a->nxloadloss, b->nxloadloss);
  TEST_EQUAL(a->startxloadlossloc, b->startxloadlossloc);
  TEST_EQUAL(a->startxloadlosslocglob, b->startxloadlosslocglob);
  TEST_EQUAL(a->nx, b->nx);
  TEST_EQUAL(a->nconeq, b->nconeq);
  TEST_EQUAL(a->nconineq, b->nconineq);
  TEST_EQUAL(a->startloc, b->startloc);
  TEST_EQUAL(a->startlocglob, b->startlocglob);
  TEST_EQUAL(a->starteqloc, b->starteqloc);
  TEST_EQUAL(a->startineqloc, b->startineqloc);

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(check_equal_gen)(PSGEN a, PSGEN b) {
  TEST_EQUAL(a->bus_i, b->bus_i);
  TEST_EQUAL(Strip(a->i), Strip(b->i));
  TEST_EQUAL(Strip(a->id), Strip(b->id));
  TEST_EQUAL(a->pg, b->pg);
  TEST_EQUAL(a->qg, b->qg);
  TEST_EQUAL(a->qt, b->qt);
  TEST_EQUAL(a->qb, b->qb);
  TEST_EQUAL(a->vs, b->vs);
  TEST_EQUAL(a->ireg, b->ireg);
  TEST_EQUAL(a->mbase, b->mbase);
  TEST_EQUAL(a->zr, b->zr);
  TEST_EQUAL(a->zx, b->zx);
  TEST_EQUAL(a->rt, b->rt);
  TEST_EQUAL(a->xt, b->xt);
  TEST_EQUAL(a->gtap, b->gtap);
  TEST_EQUAL(a->status, b->status);
  TEST_EQUAL(a->initial_status, b->initial_status);
  TEST_EQUAL(a->rmpct, b->rmpct);
  TEST_EQUAL(a->pt, b->pt);
  TEST_EQUAL(a->pb, b->pb);
  TEST_EQUAL(a->o1, b->o1);
  TEST_EQUAL(a->f1, b->f1);
  TEST_EQUAL(a->internal_i, b->internal_i);
  TEST_EQUAL(a->scale_gen, b->scale_gen);
  TEST_EQUAL(a->cost_model, b->cost_model);
  TEST_EQUAL(a->cost_startup, b->cost_startup);
  TEST_EQUAL(a->cost_shutdown, b->cost_shutdown);
  TEST_EQUAL(a->cost_ncoeffs, b->cost_ncoeffs);
  TEST_EQUAL(a->cost_gamma, b->cost_gamma);
  TEST_EQUAL(a->cost_beta, b->cost_beta);
  TEST_EQUAL(a->cost_alpha, b->cost_alpha);
  TEST_EQUAL(a->pc1, b->pc1);
  TEST_EQUAL(a->pc2, b->pc2);
  TEST_EQUAL(a->qc1min, b->qc1min);
  TEST_EQUAL(a->qc1max, b->qc1max);
  TEST_EQUAL(a->qc2min, b->qc2min);
  TEST_EQUAL(a->qc2max, b->qc2max);
  TEST_EQUAL(a->ramp_rate_min, b->ramp_rate_min);
  TEST_EQUAL(a->ramp_rate_10min, b->ramp_rate_10min);
  TEST_EQUAL(a->ramp_rate_30min, b->ramp_rate_30min);
  TEST_EQUAL(a->ramp_rate_min_mvar, b->ramp_rate_min_mvar);
  TEST_EQUAL(a->apf, b->apf);
  TEST_EQUAL(a->pgs, b->pgs);
  TEST_EQUAL(a->genfuel_type, b->genfuel_type);
  TEST_EQUAL(a->nxpow, b->nxpow);
  TEST_EQUAL(a->nxpset, b->nxpset);
  TEST_EQUAL(a->nxpdev, b->nxpdev);
  TEST_EQUAL(a->startxpowloc, b->startxpowloc);
  TEST_EQUAL(a->startxpsetloc, b->startxpsetloc);
  TEST_EQUAL(a->startxpdevloc, b->startxpdevloc);
  TEST_EQUAL(a->startxpowlocglob, b->startxpowlocglob);
  TEST_EQUAL(a->startxpsetlocglob, b->startxpsetlocglob);
  TEST_EQUAL(a->startxpdevlocglob, b->startxpdevlocglob);
  TEST_EQUAL(a->nx, b->nx);
  TEST_EQUAL(a->nconeq, b->nconeq);
  TEST_EQUAL(a->nconineq, b->nconineq);
  TEST_EQUAL(a->startloc, b->startloc);
  TEST_EQUAL(a->starteqloc, b->starteqloc);
  TEST_EQUAL(a->startineqloc, b->startineqloc);
  TEST_EQUAL(a->isrenewable, b->isrenewable);

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(check_equal_line)(PSLINE a, PSLINE b) {
  TEST_EQUAL(a->fbus, b->fbus);
  TEST_EQUAL(a->tbus, b->tbus);
  TEST_EQUAL(Strip(a->i), Strip(b->i));
  TEST_EQUAL(Strip(a->j), Strip(b->j));
  TEST_EQUAL(Strip(a->ckt), Strip(b->ckt));
  TEST_EQUAL(a->r, b->r);
  TEST_EQUAL(a->x, b->x);
  TEST_EQUAL(a->b, b->b);
  TEST_EQUAL(a->rateA, b->rateA);
  TEST_EQUAL(a->rateB, b->rateB);
  TEST_EQUAL(a->rateC, b->rateC);
  TEST_EQUAL(a->tapratio, b->tapratio);
  TEST_EQUAL(a->phaseshift, b->phaseshift);
  TEST_EQUAL(a->gi, b->gi);
  TEST_EQUAL(a->bi, b->bi);
  TEST_EQUAL(a->gj, b->gj);
  TEST_EQUAL(a->bj, b->bj);
  TEST_EQUAL(a->status, b->status);
  // TEST_EQUAL(a->met, b->met); // hard-coded in psreaddata
  TEST_EQUAL(a->length, b->length);
  TEST_EQUAL(a->o1, b->o1);
  TEST_EQUAL(a->f1, b->f1);
  TEST_EQUAL(a->sbase12, b->sbase12);
  TEST_EQUAL(a->internal_i, b->internal_i);
  TEST_EQUAL(a->internal_j, b->internal_j);
  TEST_EQUAL(a->yff[0], b->yff[0]);
  TEST_EQUAL(a->yff[1], b->yff[1]);
  TEST_EQUAL(a->yft[0], b->yft[0]);
  TEST_EQUAL(a->yft[1], b->yft[1]);
  TEST_EQUAL(a->ytf[0], b->ytf[0]);
  TEST_EQUAL(a->ytf[1], b->ytf[1]);
  TEST_EQUAL(a->ytt[0], b->ytt[0]);
  TEST_EQUAL(a->ytt[1], b->ytt[1]);
  TEST_EQUAL(a->bdc, b->bdc);
  TEST_EQUAL(a->pshift, b->pshift);
  TEST_EQUAL(a->pf, b->pf);
  TEST_EQUAL(a->qf, b->qf);
  TEST_EQUAL(a->pt, b->pt);
  TEST_EQUAL(a->qt, b->qt);
  TEST_EQUAL(a->sf, b->sf);
  TEST_EQUAL(a->st, b->st);
  TEST_EQUAL(a->reversed_ends, b->reversed_ends);
  TEST_EQUAL(a->kvlevel, b->kvlevel);
  TEST_EQUAL(a->areaf, b->areaf);
  TEST_EQUAL(a->areat, b->areat);
  TEST_EQUAL(a->zonef, b->zonef);
  TEST_EQUAL(a->zonet, b->zonet);
  TEST_EQUAL(a->connbuses[0], b->connbuses[0]);
  TEST_EQUAL(a->connbuses[1], b->connbuses[1]);

  TEST_EQUAL(a->isdcline, b->isdcline);
  TEST_EQUAL(a->pmin, b->pmin);
  TEST_EQUAL(a->pmax, b->pmax);
  TEST_EQUAL(a->Vf, b->Vf);
  TEST_EQUAL(a->Vt, b->Vt);
  TEST_EQUAL(a->qminf, b->qminf);
  TEST_EQUAL(a->qmaxf, b->qmaxf);
  TEST_EQUAL(a->qmint, b->qmint);
  TEST_EQUAL(a->qmaxt, b->qmaxt);
  TEST_EQUAL(a->loss0, b->loss0);
  TEST_EQUAL(a->loss1, b->loss1);
  TEST_EQUAL(a->mult_pmin, b->mult_pmin);
  TEST_EQUAL(a->mult_pmax, b->mult_pmax);
  TEST_EQUAL(a->mult_qminf, b->mult_qminf);
  TEST_EQUAL(a->mult_qmaxf, b->mult_qmaxf);
  TEST_EQUAL(a->mult_qmint, b->mult_qmint);
  TEST_EQUAL(a->mult_qmaxt, b->mult_qmaxt);
  TEST_EQUAL(a->mult_pf, b->mult_pf);
  TEST_EQUAL(a->mult_sf, b->mult_sf);
  TEST_EQUAL(a->mult_st, b->mult_st);
  TEST_EQUAL(a->subst_from, b->subst_from);
  TEST_EQUAL(a->subst_to, b->subst_to);
  TEST_EQUAL(a->startloc, b->startloc);
  TEST_EQUAL(a->startxslackloc, b->startxslackloc);
  TEST_EQUAL(a->startlocglob, b->startlocglob);
  TEST_EQUAL(a->nx, b->nx);
  TEST_EQUAL(a->nconeq, b->nconeq);
  TEST_EQUAL(a->nconineq, b->nconineq);
  TEST_EQUAL(a->startxdcloc, b->startxdcloc);
  TEST_EQUAL(a->starteqloc, b->starteqloc);
  TEST_EQUAL(a->startineqloc, b->startineqloc);

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(check_equal)(PS n1, PS n2) {
  TEST_EQUAL(n1->MVAbase, n2->MVAbase);
  TEST_EQUAL(n1->Nbus, n2->Nbus);
  TEST_EQUAL(n1->Ngen, n2->Ngen);
  TEST_EQUAL(n1->Nline, n2->Nline);
  TEST_EQUAL(n1->Nload, n2->Nload);
  TEST_EQUAL(n1->Ndcline, n2->Ndcline);
  TEST_EQUAL(n1->nbus, n2->nbus);
  TEST_EQUAL(n1->ngen, n2->ngen);
  TEST_EQUAL(n1->nline, n2->nline);
  TEST_EQUAL(n1->nload, n2->nload);
  TEST_EQUAL(n1->ndcline, n2->ndcline);
  TEST_EQUAL(n1->NlineON, n2->NlineON);
  TEST_EQUAL(n1->nlineON, n2->nlineON);
  TEST_EQUAL(n1->NdclineON, n2->NdclineON);
  TEST_EQUAL(n1->ndclineON, n2->ndclineON);
  TEST_EQUAL(n1->NgenON, n2->NgenON);
  TEST_EQUAL(n1->ngenON, n2->ngenON);
  TEST_EQUAL(n1->Nref, n2->Nref);
  TEST_EQUAL(n1->nref, n2->nref);

  TEST_EQUAL(n1->ngencoal, n2->ngencoal);
  TEST_EQUAL(n1->ngenwind, n2->ngenwind);
  TEST_EQUAL(n1->ngensolar, n2->ngensolar);
  TEST_EQUAL(n1->ngenng, n2->ngenng);
  TEST_EQUAL(n1->ngennuclear, n2->ngennuclear);
  TEST_EQUAL(n1->ngenhydro, n2->ngenhydro);
  TEST_EQUAL(n1->ngenundefined, n2->ngenundefined);
  TEST_EQUAL(n1->ngenrenew, n2->ngenrenew);

  TEST_EQUAL(n1->nisolated_buses, n2->nisolated_buses);
  TEST_EQUAL(n1->Nisolated_buses, n2->Nisolated_buses);
  for (int i = 0; i < n1->Nisolated_buses; ++i) {
    TEST_EQUAL(n1->isolated_buses[i], n2->isolated_buses[i]);
  }

  for (int i = 0; i < n1->Nbus; ++i) {
    TEST(check_equal_bus(&n1->bus[i], &n2->bus[i]));
  }
  for (int i = 0; i < n1->Nload; ++i) {
    TEST(check_equal_load(&n1->load[i], &n2->load[i]));
  }
  for (int i = 0; i < n1->Ngen; ++i) {
    TEST(check_equal_gen(&n1->gen[i], &n2->gen[i]));
  }
  for (int i = 0; i < n1->Nline; ++i) {
    TEST(check_equal_line(&n1->line[i], &n2->line[i]));
  }

  TEST_EQUAL(n1->refct, n2->refct);
  TEST_EQUAL(n1->maxbusnum, n2->maxbusnum);
  TEST_EQUAL(n1->ndiff, n2->ndiff);

  TEST_EQUAL(Strip(n1->net_file_name), Strip(n2->net_file_name));
  TEST_EQUAL(Strip(n1->gic_file_name), Strip(n2->gic_file_name));
  TEST_EQUAL(n1->gic_file_set, n2->gic_file_set);

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(check_bus_ids)(const exago::psse::Network &nw) {
  for (std::size_t i = 0; i < nw.buses.size(); ++i) {
    const auto &bus = nw.buses[i];
    TEST_EQUAL(nw.bus_mapping.GetInternalIndex(bus.i), i);
    TEST_EQUAL(nw.bus_mapping.GetInternalIndex(bus.name), i);
    TEST_EQUAL(nw.bus_mapping.GetBusName(bus.i), bus.name);
    TEST_EQUAL(nw.bus_mapping.GetBusNumber(bus.name), bus.i);
    TEST_EQUAL(&bus, &nw.bus_mapping.GetBus(bus.i));
    TEST_EQUAL(&bus, &nw.bus_mapping.GetBus(bus.name));
  }
  for (auto &&load : nw.loads) {
    const auto &bus = nw.bus_mapping.GetBus(load.i);
    TEST_EQUAL(&bus, &nw.bus_mapping.GetBus(load.i_bus_name));
    TEST_EQUAL(bus.i, load.i);
    TEST_EQUAL(bus.name, load.i_bus_name);
  }
  for (auto &&shunt : nw.fixed_bus_shunts) {
    const auto &bus = nw.bus_mapping.GetBus(shunt.i);
    TEST_EQUAL(&bus, &nw.bus_mapping.GetBus(shunt.i_bus_name));
    TEST_EQUAL(bus.i, shunt.i);
    TEST_EQUAL(bus.name, shunt.i_bus_name);
  }
  for (auto &&gen : nw.generators) {
    const auto &bus = nw.bus_mapping.GetBus(gen.i);
    TEST_EQUAL(&bus, &nw.bus_mapping.GetBus(gen.i_bus_name));
    TEST_EQUAL(bus.i, gen.i);
    TEST_EQUAL(bus.name, gen.i_bus_name);
    if (gen.ireg == 0) {
      continue;
    }
    const auto &regbus = nw.bus_mapping.GetBus(gen.ireg);
    TEST_EQUAL(&regbus, &nw.bus_mapping.GetBus(gen.ireg_bus_name));
    TEST_EQUAL(regbus.i, gen.ireg);
    TEST_EQUAL(regbus.name, gen.ireg_bus_name);
  }
  for (auto &&br : nw.branches) {
    const auto &ibus = nw.bus_mapping.GetBus(br.i);
    TEST_EQUAL(&ibus, &nw.bus_mapping.GetBus(br.i_bus_name));
    TEST_EQUAL(ibus.i, br.i);
    TEST_EQUAL(ibus.name, br.i_bus_name);
    const auto &jbus = nw.bus_mapping.GetBus(br.j);
    TEST_EQUAL(&jbus, &nw.bus_mapping.GetBus(br.j_bus_name));
    TEST_EQUAL(jbus.i, br.j);
    TEST_EQUAL(jbus.name, br.j_bus_name);
  }
  for (auto &&tr : nw.transformers) {
    const auto &ibus = nw.bus_mapping.GetBus(tr.i);
    TEST_EQUAL(&ibus, &nw.bus_mapping.GetBus(tr.i_bus_name));
    TEST_EQUAL(ibus.i, tr.i);
    TEST_EQUAL(ibus.name, tr.i_bus_name);
    const auto &jbus = nw.bus_mapping.GetBus(tr.j);
    TEST_EQUAL(&jbus, &nw.bus_mapping.GetBus(tr.j_bus_name));
    TEST_EQUAL(jbus.i, tr.j);
    TEST_EQUAL(jbus.name, tr.j_bus_name);
    if (tr.k == 0) {
      continue;
    }
    const auto &kbus = nw.bus_mapping.GetBus(tr.k);
    TEST_EQUAL(&kbus, &nw.bus_mapping.GetBus(tr.k_bus_name));
    TEST_EQUAL(kbus.i, tr.k);
    TEST_EQUAL(kbus.name, tr.k_bus_name);
  }

  for (auto &sh : nw.switched_shunts) {
    const auto &bus = nw.bus_mapping.GetBus(sh.i);
    TEST_EQUAL(&bus, &nw.bus_mapping.GetBus(sh.i_bus_name));
    TEST_EQUAL(bus.i, sh.i);
    TEST_EQUAL(bus.name, sh.i_bus_name);
    if (sh.swrem == 0) {
      continue;
    }
    const auto &swrembus = nw.bus_mapping.GetBus(sh.swrem);
    TEST_EQUAL(&swrembus, &nw.bus_mapping.GetBus(sh.swrem_bus_name));
    TEST_EQUAL(swrembus.i, sh.swrem);
    TEST_EQUAL(swrembus.name, sh.swrem_bus_name);
  }

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

TEST_FUNCTION(ieee9bus_v33)() {
  std::string filename{"ieee9bus_v33.raw"};

  auto psh = ReadPSData(filename);
  TEST(check_ps_data_ieee9bus(psh.ps));

  auto nw = exago::psse::ParseNetwork(filename);
  auto nw_psh = NetworkToPS(nw);
  TEST(check_ps_data_ieee9bus(nw_psh.ps));

  TEST(check_equal(psh.ps, nw_psh.ps));

  TEST(check_bus_ids(nw));

  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

// TEST_FUNCTION(check_network_ieee9bus_shunts)(const exago::psse::Network &nw)
// {
//   TEST_EQUAL(nw.case_id.ic, 0);
//   TEST_EQUAL(nw.case_id.sbase, 100.0);
//   TEST_EQUAL(nw.case_id.rev, 33);

//   const auto &buses = nw.buses;
//   TEST_EQUAL(buses[0].i, 1);
//   TEST_EQUAL(buses[0].name, "WINNSBORO 0");
//   TEST_EQUAL(buses[0].evlo, 0.9);
//   TEST_EQUAL(buses[1].i, 500);
//   TEST_EQUAL(buses[1].name, "MC CORMICK 0");
//   TEST_EQUAL(buses[1].evlo, 0.9);

//   const auto& loads = nw.loads;
//   TEST_EQUAL(loads[0].i, 2);
//   TEST_EQUAL(loads[0].scale, 1);
//   TEST_EQUAL(loads[0].intrpt, 0);
//   TEST_EQUAL(loads[1].i, 500);
//   TEST_EQUAL(loads[1].scale, 1);
//   TEST_EQUAL(loads[1].intrpt, 0);

//   const auto& fixed_shunts = nw.fixed_bus_shunts;
//   TEST_EQUAL(fixed_shunts[0].i, 320);
//   TEST_EQUAL(fixed_shunts[0].bl, 48.274);
//   TEST_EQUAL(fixed_shunts[1].i, 784);
//   TEST_EQUAL(fixed_shunts[1].bl, 19.939);

//   const auto& gens = nw.generators;
//   TEST_EQUAL(gens[0].i, 9);
//   TEST_EQUAL(gens[0].wpf, 1.0);
//   TEST_EQUAL(gens[1].i, 498);
//   TEST_EQUAL(gens[1].wpf, 1.0);

//   const auto& branches = nw.branches;
//   TEST_EQUAL(branches[0].i, 2);
//   TEST_EQUAL(branches[0].met, 1);
//   TEST_EQUAL(branches[0].owners[3].fraction, 1.0);
//   TEST_EQUAL(branches[1].i, 500);
//   TEST_EQUAL(branches[1].met, 1);
//   TEST_EQUAL(branches[1].owners[3].fraction, 1.0);

//   const auto& transformers = nw.transformers;
//   TEST_EQUAL(transformers[0].i, 8);
//   TEST_EQUAL(transformers[1].i, 190);
//   TEST_EQUAL(transformers[2].i, 498);
//   TEST(transformers[0].vecgrp.empty());
//   TEST(transformers[1].vecgrp.empty());
//   TEST(transformers[2].vecgrp.empty());
//   TEST_EQUAL(transformers[0].imp12.r, 3.55062e-4);
//   TEST_EQUAL(transformers[1].imp12.r, 7.65222e-4);
//   TEST_EQUAL(transformers[2].imp12.r, 4.00610e-3);
//   TEST_EQUAL(transformers[0].imp12.sbase, 100.0);
//   TEST_EQUAL(transformers[1].imp12.sbase, 100.0);
//   TEST_EQUAL(transformers[2].imp12.sbase, 100.0);
//   // TEST_EQUAL(transformers[0].anstar, 0.0);
//   // TEST_EQUAL(transformers[1].anstar, -65.843144);
//   // TEST_EQUAL(transformers[2].anstar, 0.0);
//   TEST_EQUAL(transformers[0].windings[0].windv, 1.0);
//   TEST_EQUAL(transformers[1].windings[0].windv, 1.0);
//   TEST_EQUAL(transformers[2].windings[0].windv, 1.0);
//   TEST_EQUAL(transformers[0].windings[0].cnxa, 0.0);
//   TEST_EQUAL(transformers[1].windings[0].cnxa, 0.0);
//   TEST_EQUAL(transformers[2].windings[0].cnxa, 0.0);
//   TEST_EQUAL(transformers[0].windings[1].windv, 1.0);
//   TEST_EQUAL(transformers[1].windings[1].windv, 1.0);
//   TEST_EQUAL(transformers[2].windings[1].windv, 1.0);
//   TEST_EQUAL(transformers[0].windings[1].nomv, 345.0);
//   TEST_EQUAL(transformers[1].windings[1].nomv, 115.0);
//   TEST_EQUAL(transformers[2].windings[1].nomv, 138.0);
//   TEST_EQUAL(transformers[0].windings[1].cnxa, 0.0);
//   TEST_EQUAL(transformers[1].windings[1].cnxa, 0.0);
//   TEST_EQUAL(transformers[2].windings[1].cnxa, 0.0);
//   TEST_EQUAL(transformers[0].windings[2].windv, 1.0);
//   TEST_EQUAL(transformers[1].windings[2].windv, 1.0);
//   TEST_EQUAL(transformers[2].windings[2].windv, 1.0);
//   TEST_EQUAL(transformers[0].windings[2].cnxa, 0.0);
//   TEST_EQUAL(transformers[1].windings[2].cnxa, 0.0);
//   TEST_EQUAL(transformers[2].windings[2].cnxa, 0.0);
//   TEST_EQUAL(transformers[0].ckt, "1");

//   // const auto& areaInterchanges = nw.areaInterchanges;
//   // TEST_EQUAL(areaInterchanges[0].i, 1);
//   // TEST_EQUAL(areaInterchanges[0].arname, "SouthCarolin");

//   // twoTerminalDC = nw.twoTerminalDCLines;
//   // TEST_EQUAL(twoTerminalDC[0].name,"DC Line 1");
//   // TEST_EQUAL(twoTerminalDC[1].name,"DC Line 1");
//   // TEST_EQUAL(twoTerminalDC[0].cccacc, 0.0);
//   // TEST_EQUAL(twoTerminalDC[1].cccacc, 0.0);
//   // TEST_EQUAL(twoTerminalDC[0].ipr, 2060653);
//   // TEST_EQUAL(twoTerminalDC[1].ipr, 3008030);
//   // TEST_EQUAL(twoTerminalDC[0].xcapr, 0.0);
//   // TEST_EQUAL(twoTerminalDC[1].xcapr, 0.0);
//   // TEST_EQUAL(twoTerminalDC[0].ipi, 66353);
//   // TEST_EQUAL(twoTerminalDC[1].ipi, 61477);
//   // TEST_EQUAL(twoTerminalDC[0].xcapi, 0.0);
//   // TEST_EQUAL(twoTerminalDC[1].xcapi, 0.0);

//   // TEST(nw.vscDCLines.empty());

//   // const auto& impedanceCorrections = nw.impedanceCorrections;
//   // TEST_EQUAL(impedanceCorrections[0].i, 1);
//   // TEST_EQUAL(impedanceCorrections[1].i, 2);
//   // TEST_EQUAL(impedanceCorrections[2].i, 3);
//   // TEST_EQUAL(impedanceCorrections[0].f6, 1.03);
//   // TEST_EQUAL(impedanceCorrections[1].f6, 0.0);
//   // TEST_EQUAL(impedanceCorrections[2].f6, 1.41);
//   // TEST_EQUAL(impedanceCorrections[0].f11, 0.0);
//   // TEST_EQUAL(impedanceCorrections[1].f11, 0.0);
//   // TEST_EQUAL(impedanceCorrections[2].f11, 0.0);

//   // TODO: other elements

//   const auto& switched_shunts = nw.switched_shunts;
//   TEST_EQUAL(switched_shunts[0].i, 27);
//   TEST_EQUAL(switched_shunts[1].i, 491);
//   TEST_EQUAL(switched_shunts[0].adjm, 0);
//   TEST_EQUAL(switched_shunts[1].adjm, 0);
//   TEST_EQUAL(switched_shunts[0].stat, 1);
//   TEST_EQUAL(switched_shunts[1].stat, 1);
//   TEST_EQUAL(switched_shunts[0].blocks[0].n, 1);
//   TEST_EQUAL(switched_shunts[1].blocks[0].n, 1);
//   TEST_EQUAL(switched_shunts[0].blocks[0].b, 80.0);
//   TEST_EQUAL(switched_shunts[1].blocks[0].b, 50.0);
//   TEST_EQUAL(switched_shunts[0].blocks[7].n, 0);
//   TEST_EQUAL(switched_shunts[1].blocks[7].n, 0);
//   TEST_EQUAL(switched_shunts[0].blocks[7].b, 0.0);
//   TEST_EQUAL(switched_shunts[1].blocks[7].b, 0.0);

//   TEST_FUNCTION_RETURN;
// }
// END_TEST_FUNCTION

// TEST_FUNCTION(ieee9bus_v34_shunts)() {
//   std::string filename{"ieee9bus_v34_shunts.raw"};

//   auto nw = exago::psse::ParseNetwork(filename);
//   TEST(check_network_ieee9bus_shunts(nw));
//   TEST(check_bus_ids(nw));

//   TEST_FUNCTION_RETURN;
// }
// END_TEST_FUNCTION

TEST_FUNCTION(driver)() {
  TEST(ieee9bus_v33());
  // TEST(ieee9bus_v34_shunts());
  TEST_FUNCTION_RETURN;
}
END_TEST_FUNCTION

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  driver();

  MPI_Finalize();
  return TEST_RESULT;
}
