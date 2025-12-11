#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#include <psse.hpp>
#include <psimpl.h>

namespace exago {
namespace psse {

double deg2rad(double deg) { return deg * M_PI / 180.0; }

std::string strip(std::string str) {
  auto notspace = [](char c) { return !std::isspace(c); };
  str.erase(begin(str), std::find_if(begin(str), end(str), notspace));
  str.erase(std::find_if(rbegin(str), rend(str), notspace).base(), end(str));
  return str;
}

std::string read_line(std::istream &is) {
  std::string line;
  std::getline(is, line);
  return strip(line);
}

class SQString {
public:
  operator std::string() const { return s_; }

private:
  friend std::istream &operator>>(std::istream &is, SQString &sqs) {
    // Remove leading whitespace
    while (std::isspace(is.peek())) {
      is.get();
    }
    // Check for and remove opening quote
    if (is.peek() != '\'') {
      throw std::runtime_error("First character expected to be single quote");
    }
    is.get();
    // Get line to closing quote
    std::getline(is, sqs.s_, '\'');

    return is;
  }

  std::string s_;
};

class LineItemStream : public std::istream {
public:
  LineItemStream() = delete;
  LineItemStream(std::istream &is) : is_(&is) { next_line(); }

  operator bool() const { return static_cast<bool>(ss_); }

  std::size_t size() const noexcept { return size_; }

  std::string str() const { return ss_.str(); }

  bool starts_with(const std::string &sub) const {
    return ss_.str().find(sub) == 0;
  }

  LineItemStream &next_line() {
    auto line = read_line(*is_);
    std::istringstream iss(line);
    ss_.str("");
    for (std::string item; std::getline(iss, item, ',');) {
      ss_ << item << ' ';
      ++size_;
    }
    return *this;
  }

  template <typename T> LineItemStream &operator>>(T &item) {
    ss_ >> item;
    --size_;
    return *this;
  }

private:
  std::istream *is_{nullptr};
  std::stringstream ss_;
  std::size_t size_{0};
};

CaseID parse_case_id(std::istream &is) {
  CaseID cid;
  LineItemStream lis(is);
  lis >> cid.ic >> cid.sbase >> cid.rev >> cid.xfrrat >> cid.nxfrat >>
      cid.basfrq;
  cid.extra[0] = read_line(lis);
  cid.extra[1] = read_line(is);
  cid.extra[2] = read_line(is);
  return cid;
}

void parse_record(LineItemStream &lis, Bus &bus) {
  SQString name;
  lis >> bus.i >> name >> bus.baskv >> bus.ide >> bus.area >> bus.zone >>
      bus.owner >> bus.vm >> bus.va;
  bus.name = strip(name);

  // TODO: parse remaining bus items if present
}

void parse_record(LineItemStream &lis, Load &ld) {
  SQString id;
  lis >> ld.i >> id >> ld.status >> ld.area >> ld.zone >> ld.pl >> ld.ql >>
      ld.ip >> ld.iq >> ld.yp >> ld.yq >> ld.owner;
  ld.id = strip(id);

  // TODO: parse remaining load items if present
}

void parse_record(LineItemStream &lis, FixedBusShunt sh) {
  SQString id;
  lis >> sh.i >> id >> sh.status >> sh.gl >> sh.bl;
  sh.id = strip(id);
}

void parse_record(LineItemStream &lis, Generator &gen) {
  SQString id;
  lis >> gen.i >> id >> gen.pg >> gen.qg >> gen.qt >> gen.qb >> gen.vs >>
      gen.ireg >> gen.mbase >> gen.zr >> gen.zx >> gen.rt >> gen.xt >>
      gen.gtap >> gen.stat >> gen.rmpct >> gen.pt >> gen.pb >>
      gen.owners[0].owner >> gen.owners[0].fraction;
  gen.id = strip(id);
}

void parse_record(LineItemStream &lis, Branch &br) {
  SQString ckt;
  lis >> br.i >> br.j >> ckt >> br.r >> br.x >> br.b >> br.ratea >> br.rateb >>
      br.ratec >> br.gi >> br.bi >> br.gj >> br.bj >> br.st >> br.met >>
      br.len >> br.owners[0].owner >> br.owners[0].fraction;
  br.ckt = strip(ckt);
}

Winding parse_transformer_winding(LineItemStream &lis) {
  Winding w;
  lis >> w.windv >> w.nomv >> w.ang >> w.rata >> w.ratb >> w.ratc >> w.cod >>
      w.cont >> w.rma >> w.rmi >> w.vma >> w.vmi >> w.ntp >> w.tab >> w.cr >>
      w.cx >> w.cnxa;
  return w;
}

void parse_record(LineItemStream &lis, Transformer &tr) {
  SQString ckt;
  SQString name;
  lis >> tr.i >> tr.j >> tr.k >> ckt >> tr.cw >> tr.cz >> tr.cm >> tr.mag1 >>
      tr.mag2 >> tr.nmetr >> name >> tr.stat >> tr.owners[0].owner >>
      tr.owners[0].fraction;
  tr.ckt = strip(ckt);
  tr.name = strip(name);
  if (tr.k == 0) {
    // two-winding (3 more rows)
    lis.next_line() >> tr.imp12.r >> tr.imp12.x >> tr.imp12.sbase;
    tr.windings[0] = parse_transformer_winding(lis.next_line());
    lis.next_line() >> tr.windings[1].windv >> tr.windings[1].nomv;
  } else {
    // three-winding (4 more rows)
    lis.next_line() >> tr.imp12.r >> tr.imp12.x >> tr.imp12.sbase >>
        tr.imp23.r >> tr.imp23.x >> tr.imp23.sbase >> tr.imp31.r >>
        tr.imp31.x >> tr.imp31.sbase >> tr.vmstar >> tr.anstar;
    tr.windings[0] = parse_transformer_winding(lis.next_line());
    tr.windings[1] = parse_transformer_winding(lis.next_line());
    tr.windings[2] = parse_transformer_winding(lis.next_line());
  }
}

template <typename T> std::vector<T> parse_records(std::istream &is) {
  std::vector<T> recs;
  while (is) {
    LineItemStream lis(is);
    if (lis.starts_with("0 /")) {
      break;
    }
    auto &rec = recs.emplace_back();
    parse_record(lis, rec);
  }
  return recs;
}

std::unordered_map<int, int> process_bus_ids(const std::vector<Bus>& buses) {
  std::unordered_map<int, int> id_map;
  for (std::size_t i = 0; i < buses.size(); ++i) {
    id_map.emplace(buses[i].i, i);
  }
  return id_map;
}

Network parse_network(std::istream &is) {
  auto case_id = parse_case_id(is);
  auto buses = parse_records<Bus>(is);
  auto bus_id_map = process_bus_ids(buses);
  auto loads = parse_records<Load>(is);
  auto shunts = parse_records<FixedBusShunt>(is);
  auto generators = parse_records<Generator>(is);
  auto branches = parse_records<Branch>(is);
  auto transformers = parse_records<Transformer>(is);

  Network nw{"",
             std::move(case_id),
             std::move(buses),
             std::move(bus_id_map),
             std::move(loads),
             std::move(shunts),
             std::move(generators),
             std::move(branches),
             std::move(transformers)};
  return nw;
}

Network parse_network(const std::string &filename) {
  std::ifstream is(filename);
  auto nw = parse_network(is);
  nw.file_name = filename;
  return nw;
}

PetscErrorCode convert_to_ps(PS ps, const Network &nw) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ps->NgenON = 0;
  ps->NlineON = 0;
  ps->nlineON = 0;

  strcpy(ps->net_file_name, nw.file_name.c_str());
  ps->MVAbase = nw.case_id.sbase;
  ps->maxbusnum =
      std::max_element(begin(nw.buses), end(nw.buses),
                       [](auto &&b1, auto &&b2) { return b1.i < b2.i; })
          ->i;

  // buses
  ps->Nbus = ps->nbus = nw.buses.size();
  ierr = PetscCalloc1(ps->Nbus, &ps->bus);
  CHKERRQ(ierr);
  for (int i = 0; i < ps->Nbus; ++i) {
    auto &dbus = ps->bus[i];
    const auto &sbus = nw.buses[i];
    dbus.bus_i = sbus.i;
    strcpy(dbus.name, sbus.name.c_str());
    dbus.basekV = sbus.baskv;
    dbus.ide = sbus.ide;
    dbus.area = sbus.area;
    dbus.zone = sbus.zone;
    dbus.owner = sbus.owner;
    dbus.vm = sbus.vm;
    dbus.va = deg2rad(sbus.va);
    dbus.nvhi = sbus.nvhi;
    dbus.nvlo = sbus.nvlo;
    dbus.evhi = sbus.evhi;
    dbus.evlo = sbus.evlo;

    if (dbus.ide == REF_BUS) {
      ps->Nref++;
    }
    dbus.internal_i = i;
    dbus.nload = 0;
    dbus.ngen = 0;
    dbus.ngenON = 0;
    dbus.nshunt = 0;
    dbus.Vmin = 1.1;
    dbus.Vmax = 0.9;
    dbus.gl = 0;
    dbus.bl = 0;
    dbus.qrange = 0.0;
    dbus.qmintot = 0.0;
    dbus.Pgtot = 0.0;
    dbus.MVAbasetot = 0.0;
  }

  // loads
  ps->Nload = ps->nload = nw.loads.size();
  ierr = PetscCalloc1(ps->Nload, &ps->load);
  CHKERRQ(ierr);
  for (int i = 0; i < ps->Nload; ++i) {
    auto &dload = ps->load[i];
    const auto &sload = nw.loads[i];
    dload.bus_i = sload.i;
    strcpy(dload.id, sload.id.c_str());
    dload.status = sload.status;
    dload.area = sload.area;
    dload.zone = sload.zone;
    dload.pl = sload.pl / ps->MVAbase;
    dload.ql = sload.ql / ps->MVAbase;
    dload.ip = sload.ip / ps->MVAbase;
    dload.iq = sload.iq / ps->MVAbase;
    dload.yp = sload.yp / ps->MVAbase;
    dload.yq = sload.yq / ps->MVAbase;
    dload.owner = sload.owner;
    dload.scale = sload.scale;
    dload.intrpt = sload.intrpt;

    auto bus_i = nw.bus_id_map.at(dload.bus_i);
    dload.internal_i = bus_i;
    auto &bus = ps->bus[bus_i];
    bus.lidx[bus.nload] = i;
    bus.nload++;
  }

  // shunts
  for (auto &shunt : nw.shunts) {
    if (shunt.status == 0) {
      continue;
    }
    auto bus_i = nw.bus_id_map.at(shunt.i);
    if (ps->bus[bus_i].nshunt > 0) {
      throw std::runtime_error(
          "Bus " + std::to_string(shunt.i) +
          ": more than one fixed shunt at bus not supported");
    }
    ps->bus[bus_i].nshunt++;
    ps->bus[bus_i].gl = shunt.gl / ps->MVAbase;
    ps->bus[bus_i].bl = shunt.bl / ps->MVAbase;
  }

  // generators
  ps->Ngen = ps->ngen = nw.generators.size();
  ierr = PetscCalloc1(ps->Ngen, &ps->gen);
  CHKERRQ(ierr);
  for (int i = 0; i < ps->Ngen; ++i) {
    auto &dgen = ps->gen[i];
    auto &sgen = nw.generators[i];
    dgen.bus_i = sgen.i;
    strcpy(dgen.id, sgen.id.c_str());
    dgen.pg = sgen.pg / ps->MVAbase;
    dgen.qg = sgen.qg / ps->MVAbase;
    dgen.qt = sgen.qt / ps->MVAbase;
    dgen.qb = sgen.qb / ps->MVAbase;
    dgen.vs = sgen.vs;
    dgen.ireg = sgen.ireg;
    dgen.mbase = sgen.mbase;
    dgen.zr = sgen.zr;
    dgen.zx = sgen.zx;
    dgen.rt = sgen.rt;
    dgen.xt = sgen.xt;
    dgen.gtap = sgen.gtap;
    dgen.status = sgen.stat;
    dgen.rmpct = sgen.rmpct;
    dgen.pt = sgen.pt / ps->MVAbase;
    dgen.pb = sgen.pb / ps->MVAbase;
    dgen.o1 = sgen.owners[0].owner;
    dgen.f1 = sgen.owners[0].fraction;

    dgen.initial_status = dgen.status;
    auto bus_i = nw.bus_id_map.at(dgen.bus_i);
    dgen.internal_i = bus_i;
    auto &bus = ps->bus[bus_i];
    bus.gidx[bus.ngen] = i;
    bus.ngen++;
    if (dgen.status == 1) {
      if ((dgen.pb <= dgen.pg) && (dgen.pg <= dgen.pt)) {
        dgen.pgs = dgen.pg;
      } else {
        dgen.pgs = (dgen.pb + dgen.pt) / 2.0;
      }
      bus.qrange += (dgen.qt - dgen.qb);
      bus.qmintot += dgen.qb;
      bus.Pgtot += PetscAbsScalar(dgen.pg);
      bus.MVAbasetot += dgen.mbase;
      if (dgen.vs != bus.vm) {
        throw std::runtime_error(
            "Generator " + sgen.id +
            ": set point voltage different from bus voltage magnitude");
      }
      bus.ngenON++;
      ps->NgenON++;
    }
  }

  // lines
  ps->Nline = ps->nline = nw.branches.size() + nw.transformers.size();
  ierr = PetscCalloc1(ps->Nline, &ps->line);
  CHKERRQ(ierr);
  auto configure_line = [&nw](auto &line) {
    line.internal_i = nw.bus_id_map.at(line.fbus);
    line.internal_j = nw.bus_id_map.at(line.tbus);
    line.tapratio = 1.0;
    line.phaseshift = 0.0;

    auto R = line.r;
    auto X = line.x;
    auto Bc = line.b;

    auto Zm = R * R + X * X;
    auto G = R / Zm;
    auto B = -X / Zm;

    auto tap = line.tapratio;
    auto shift = line.phaseshift;
    auto tap2 = tap * tap;
    auto tapr = tap * cos(shift);
    auto tapi = tap * sin(shift);

    line.yff[0] = G / tap2;
    line.yff[1] = (B + Bc / 2.0) / tap2;

    line.yft[0] = -(G * tapr - B * tapi) / tap2;
    line.yft[1] = -(B * tapr + G * tapi) / tap2;

    line.ytf[0] = -(G * tapr + B * tapi) / tap2;
    line.ytf[1] = -(B * tapr - G * tapi) / tap2;

    line.ytt[0] = G;
    line.ytt[1] = B + Bc / 2.0;
  };
  int nbranch = nw.branches.size();
  int ntrans = nw.transformers.size();
  for (int i = 0; i < nbranch; ++i) {
    auto &dline = ps->line[i];
    auto &sline = nw.branches[i];

    dline.fbus = sline.i;
    dline.tbus = sline.j;
    strcpy(dline.ckt, sline.ckt.c_str());
    dline.r = sline.r;
    dline.x = sline.x;
    dline.b = sline.b;
    dline.rateA = sline.ratea;
    dline.rateB = sline.rateb;
    dline.rateC = sline.ratec;
    dline.gi = sline.gi;
    dline.bi = sline.bi;
    dline.gj = sline.gj;
    dline.bj = sline.bj;
    dline.status = sline.st;
    dline.met = sline.met;
    dline.length = sline.len;
    dline.o1 = sline.owners[0].owner;
    dline.f1 = sline.owners[0].fraction;

    (dline.rateA == 0.0) && (dline.rateA = PETSC_INFINITY);
    if (dline.status == 1) {
      ps->NlineON++;
    }

    configure_line(dline);

    dline.subst_from = dline.subst_to = nullptr;
  }
  for (int i = 0; i < ntrans; ++i) {
    auto &dline = ps->line[i + nbranch];
    auto &sline = nw.transformers[i];

    dline.fbus = sline.i;
    dline.tbus = sline.j;
    // k skipped
    strcpy(dline.ckt, sline.ckt.c_str());
    // cw, cz, cm, mag1, mag2, nmetr, name skipped
    dline.status = sline.stat;
    dline.o1 = sline.owners[0].owner;
    dline.f1 = sline.owners[0].fraction;

    dline.r = sline.imp12.r;
    dline.x = sline.imp12.x;
    dline.sbase12 = sline.imp12.sbase;

    dline.tapratio = sline.windings[0].windv;
    // nomv1 skipped
    dline.phaseshift = sline.windings[0].ang;
    dline.rateA = sline.windings[0].rata;
    dline.rateB = sline.windings[0].ratb;
    dline.rateC = sline.windings[0].ratc;
    // the rest are skipped

    configure_line(dline);
  }

  PetscFunctionReturn(0);
}

} // namespace psse
} // namespace exago
