#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#include <psse.hpp>
#include <psimpl.h>

namespace exago {
namespace psse {

std::string Strip(std::string str) {
  auto notspace = [](char c) { return !std::isspace(c); };
  str.erase(begin(str), std::find_if(begin(str), end(str), notspace));
  str.erase(std::find_if(rbegin(str), rend(str), notspace).base(), end(str));
  return str;
}

std::string ReadLine(std::istream &is) {
  std::string line;
  std::getline(is, line);
  return Strip(line);
}

std::istream &SkipLeadingWhitespace(std::istream &is) {
  while (std::isspace(is.peek())) {
    is.get();
  }
  return is;
}

class QuoteStringParse {
public:
  operator std::string() const { return s_; }

private:
  friend std::istream &operator>>(std::istream &is, QuoteStringParse &qs) {
    SkipLeadingWhitespace(is);

    // Check for and remove opening quote
    char qc = is.peek();
    if (!(qc == '\'' || qc == '\"')) {
      throw std::runtime_error(
          "First character expected to be single or double quote");
    }
    is.get();

    // Get line to closing quote
    std::getline(is, qs.s_, qc);

    return is;
  }

  std::string s_;
};

class IntOrStringParse {
public:
  int GetInt() const { return i_; }
  const std::string &GetString() const { return s_; }

private:
  friend std::istream &operator>>(std::istream &is, IntOrStringParse &p) {
    SkipLeadingWhitespace(is);

    // Check for opening quote
    auto qc = is.peek();
    if (qc == '\'' || qc == '\"') {
      QuoteStringParse qs;
      is >> qs;
      p.s_ = Strip(qs);
    } else {
      is >> p.i_;
    }

    return is;
  }

  int i_{0};
  std::string s_;
};

class LineItemStream : public std::istream {
public:
  LineItemStream() = delete;
  LineItemStream(std::istream &is) : is_(&is) { NextLine(); }

  operator bool() const { return static_cast<bool>(ss_); }

  std::size_t Size() const noexcept { return size_; }

  bool StartsWith(const std::string &sub) const {
    return ss_.str().find(sub) == 0;
  }

  LineItemStream &NextLine() {
    auto line = ReadLine(*is_);
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

CaseID ParseCaseID(std::istream &is) {
  CaseID cid;
  LineItemStream lis(is);
  lis >> cid.ic >> cid.sbase >> cid.rev >> cid.xfrrat >> cid.nxfrat >>
      cid.basfrq;
  cid.extra[0] = ReadLine(lis);
  cid.extra[1] = ReadLine(is);
  cid.extra[2] = ReadLine(is);
  return cid;
}

void ParseRecord(LineItemStream &lis, Bus &bus) {
  QuoteStringParse name;
  lis >> bus.i >> name >> bus.baskv >> bus.ide >> bus.area >> bus.zone >>
      bus.owner >> bus.vm >> bus.va;
  bus.name = Strip(name);
  if (bus.name.empty()) {
    bus.name = "BUS" + std::to_string(bus.i);
  }
  // TODO: parse remaining bus items if present
}

void ParseRecord(LineItemStream &lis, Load &ld) {
  IntOrStringParse i;
  QuoteStringParse id;
  lis >> i >> id >> ld.status >> ld.area >> ld.zone >> ld.pl >> ld.ql >>
      ld.ip >> ld.iq >> ld.yp >> ld.yq >> ld.owner;
  ld.i = i.GetInt();
  ld.i_bus_name = i.GetString();
  ld.id = Strip(id);

  // TODO: parse remaining load items if present
}

void ParseRecord(LineItemStream &lis, FixedBusShunt sh) {
  IntOrStringParse i;
  QuoteStringParse id;
  lis >> i >> id >> sh.status >> sh.gl >> sh.bl;
  sh.i = i.GetInt();
  sh.i_bus_name = i.GetString();
  sh.id = Strip(id);
}

void ParseRecord(LineItemStream &lis, Generator &gen) {
  IntOrStringParse i;
  IntOrStringParse ireg;
  QuoteStringParse id;
  lis >> i >> id >> gen.pg >> gen.qg >> gen.qt >> gen.qb >> gen.vs >> ireg >>
      gen.mbase >> gen.zr >> gen.zx >> gen.rt >> gen.xt >> gen.gtap >>
      gen.stat >> gen.rmpct >> gen.pt >> gen.pb >> gen.owners[0].owner >>
      gen.owners[0].fraction;
  gen.i = i.GetInt();
  gen.i_bus_name = i.GetString();
  gen.ireg = ireg.GetInt();
  gen.ireg_bus_name = ireg.GetString();
  gen.id = Strip(id);
}

void ParseRecord(LineItemStream &lis, Branch &br) {
  IntOrStringParse i;
  IntOrStringParse j;
  QuoteStringParse ckt;
  lis >> i >> j >> ckt >> br.r >> br.x >> br.b >> br.ratea >> br.rateb >>
      br.ratec >> br.gi >> br.bi >> br.gj >> br.bj >> br.st >> br.met >>
      br.len >> br.owners[0].owner >> br.owners[0].fraction;
  br.i = i.GetInt();
  br.i_bus_name = i.GetString();
  br.j = j.GetInt();
  br.j_bus_name = j.GetString();
  br.ckt = Strip(ckt);
}

Winding parse_transformer_winding(LineItemStream &lis) {
  Winding w;
  lis >> w.windv >> w.nomv >> w.ang >> w.rata >> w.ratb >> w.ratc >> w.cod >>
      w.cont >> w.rma >> w.rmi >> w.vma >> w.vmi >> w.ntp >> w.tab >> w.cr >>
      w.cx >> w.cnxa;
  return w;
}

void ParseRecord(LineItemStream &lis, Transformer &tr) {
  IntOrStringParse i;
  IntOrStringParse j;
  IntOrStringParse k;
  QuoteStringParse ckt;
  QuoteStringParse name;
  lis >> i >> j >> k >> ckt >> tr.cw >> tr.cz >> tr.cm >> tr.mag1 >> tr.mag2 >>
      tr.nmetr >> name >> tr.stat >> tr.owners[0].owner >>
      tr.owners[0].fraction;
  tr.i = i.GetInt();
  tr.i_bus_name = i.GetString();
  tr.j = j.GetInt();
  tr.j_bus_name = j.GetString();
  tr.k = k.GetInt();
  tr.k_bus_name = k.GetString();
  tr.ckt = Strip(ckt);
  tr.name = Strip(name);
  if (tr.k == 0) {
    // two-winding (3 more rows)
    lis.NextLine() >> tr.imp12.r >> tr.imp12.x >> tr.imp12.sbase;
    tr.windings[0] = parse_transformer_winding(lis.NextLine());
    lis.NextLine() >> tr.windings[1].windv >> tr.windings[1].nomv;
  } else {
    // three-winding (4 more rows)
    lis.NextLine() >> tr.imp12.r >> tr.imp12.x >> tr.imp12.sbase >>
        tr.imp23.r >> tr.imp23.x >> tr.imp23.sbase >> tr.imp31.r >>
        tr.imp31.x >> tr.imp31.sbase >> tr.vmstar >> tr.anstar;
    tr.windings[0] = parse_transformer_winding(lis.NextLine());
    tr.windings[1] = parse_transformer_winding(lis.NextLine());
    tr.windings[2] = parse_transformer_winding(lis.NextLine());
  }
}

void ParseRecord(LineItemStream &, AreaInterchange &) {}

void ParseRecord(LineItemStream &, TwoTerminalDCLine &) {}

void ParseRecord(LineItemStream &, VSCDCLine &) {}

void ParseRecord(LineItemStream &, ImpedanceCorrection &) {}

void ParseRecord(LineItemStream &, MultiTerminalDCLine &) {}

void ParseRecord(LineItemStream &, MultiSectionLineGroup &) {}

void ParseRecord(LineItemStream &, Zone &) {}

void ParseRecord(LineItemStream &, InterAreaTransfer &) {}

void ParseRecord(LineItemStream &, Owner &) {}

void ParseRecord(LineItemStream &, FACTSDevice &) {}

void ParseRecord(LineItemStream &lis, SwitchedShunt &sh) {
  IntOrStringParse i;
  IntOrStringParse swrem;
  QuoteStringParse rmidnt;
  lis >> i >> sh.modsw >> sh.adjm >> sh.stat >> sh.vswhi >> sh.vswlo >> swrem >>
      sh.rmpct >> rmidnt >> sh.binit >> sh.blocks[0].n >> sh.blocks[0].b;
  sh.i = i.GetInt();
  sh.i_bus_name = i.GetString();
  sh.swrem = swrem.GetInt();
  sh.swrem_bus_name = swrem.GetString();
  sh.rmidnt = Strip(rmidnt);
}

void ParseRecord(LineItemStream &, GNEDevice &) {}

template <typename T> std::vector<T> ParseRecords(std::istream &is) {
  std::vector<T> recs;
  while (is && is.peek() != 'Q') {
    LineItemStream lis(is);
    if (lis.StartsWith("0 /")) {
      break;
    }
    if (lis.StartsWith("@!")) {
      continue;
    }
    auto &rec = recs.emplace_back();
    ParseRecord(lis, rec);
  }
  return recs;
}

BusMapping::BusMapping(std::vector<Bus> &buses) : buses_(buses) {
  for (std::size_t i = 0; i < buses_.size(); ++i) {
    {
      auto [_, ins] = id_map_.emplace(buses_[i].i, i);
      if (!ins) {
        throw std::runtime_error("Bus numbers non-unique");
      }
    }
    {
      auto [_, ins] = name_map_.emplace(buses_[i].name, i);
      if (!ins) {
        throw std::runtime_error("Bus names non-unique");
      }
    }
  }
}

bool BusMapping::HasBus(std::size_t bus_number) const {
  return (id_map_.count(bus_number) > 0);
}

bool BusMapping::HasBus(const std::string &bus_name) const {
  return (name_map_.count(bus_name) > 0);
}

std::size_t BusMapping::GetInternalIndex(std::size_t bus_number) const {
  return id_map_.at(bus_number);
}

std::size_t BusMapping::GetInternalIndex(const std::string &bus_name) const {
  return name_map_.at(bus_name);
}

std::size_t BusMapping::GetBusNumber(const std::string &bus_name) const {
  return buses_[GetInternalIndex(bus_name)].i;
}

const std::string &BusMapping::GetBusName(std::size_t bus_number) const {
  return buses_[GetInternalIndex(bus_number)].name;
}

const Bus &BusMapping::GetBus(const std::string &bus_name) const {
  return buses_[GetInternalIndex(bus_name)];
}

const Bus &BusMapping::GetBus(std::size_t bus_number) const {
  return buses_[GetInternalIndex(bus_number)];
}

void BusMapping::Resolve(std::size_t &bus_number, std::string &bus_name,
                         BusMapping::Optional optional) const {
  if (bus_number == 0 && bus_name.empty()) {
    if (optional) {
      return;
    } else {
      throw std::runtime_error("Cannot resolve bus id");
    }
  }
  if (bus_number == 0) {
    if (!HasBus(bus_name)) {
      throw std::runtime_error("Bus \'" + bus_name + "\' does not exist");
    }
    bus_number = GetBusNumber(bus_name);
  }
  if (bus_name.empty()) {
    if (!HasBus(bus_number)) {
      throw std::runtime_error("Bus " + std::to_string(bus_number) +
                               " does not exist");
    }
    bus_name = GetBusName(bus_number);
  }
}

void Network::ResolveBusIds() {
  for (auto &load : loads) {
    bus_mapping.Resolve(load.i, load.i_bus_name);
  }
  for (auto &shunt : fixed_bus_shunts) {
    bus_mapping.Resolve(shunt.i, shunt.i_bus_name);
  }
  for (auto &gen : generators) {
    bus_mapping.Resolve(gen.i, gen.i_bus_name);
    bus_mapping.Resolve(gen.ireg, gen.ireg_bus_name,
                        BusMapping::Optional{true});
  }
  for (auto &br : branches) {
    bus_mapping.Resolve(br.i, br.i_bus_name);
    bus_mapping.Resolve(br.j, br.j_bus_name);
  }
  for (auto &tr : transformers) {
    bus_mapping.Resolve(tr.i, tr.i_bus_name);
    bus_mapping.Resolve(tr.j, tr.j_bus_name);
    bus_mapping.Resolve(tr.k, tr.k_bus_name, BusMapping::Optional{true});
  }

  for (auto &sh : switched_shunts) {
    bus_mapping.Resolve(sh.i, sh.i_bus_name);
    bus_mapping.Resolve(sh.swrem, sh.swrem_bus_name);
  }
}

Network::Network(CaseID &&cid, std::vector<Bus> &&bus, std::vector<Load> &&load,
                 std::vector<FixedBusShunt> &&fbshunt,
                 std::vector<Generator> &&gen, std::vector<Branch> &&branch,
                 std::vector<Transformer> &&trans,
                 std::vector<SwitchedShunt> &&swshunt)
    : case_id(std::move(cid)), buses(std::move(bus)), bus_mapping(buses),
      loads(std::move(load)), fixed_bus_shunts(std::move(fbshunt)),
      generators(std::move(gen)), branches(std::move(branch)),
      transformers(std::move(trans)), switched_shunts(std::move(swshunt)) {
  ResolveBusIds();
}

Network ParseNetwork(std::istream &is) {
  auto case_id = ParseCaseID(is);
  auto buses = ParseRecords<Bus>(is);
  auto loads = ParseRecords<Load>(is);
  auto fixed_bus_shunts = ParseRecords<FixedBusShunt>(is);
  auto generators = ParseRecords<Generator>(is);
  auto branches = ParseRecords<Branch>(is);
  auto transformers = ParseRecords<Transformer>(is);
  // { TODO
  auto area_interchanges = ParseRecords<AreaInterchange>(is);
  auto two_terminal_dc_lines = ParseRecords<TwoTerminalDCLine>(is);
  auto vsc_dc_lines = ParseRecords<VSCDCLine>(is);
  auto impedance_corrections = ParseRecords<ImpedanceCorrection>(is);
  auto multi_terminal_dc_lines = ParseRecords<MultiTerminalDCLine>(is);
  auto multi_section_line_groups = ParseRecords<MultiSectionLineGroup>(is);
  auto zones = ParseRecords<Zone>(is);
  auto inter_area_transfers = ParseRecords<InterAreaTransfer>(is);
  auto owners = ParseRecords<Owner>(is);
  auto facts_devices = ParseRecords<FACTSDevice>(is);
  // }
  auto switched_shunts = ParseRecords<SwitchedShunt>(is);

  // { TODO
  auto gne_devices = ParseRecords<GNEDevice>(is);
  // }

  Network nw{std::move(case_id),      std::move(buses),
             std::move(loads),        std::move(fixed_bus_shunts),
             std::move(generators),   std::move(branches),
             std::move(transformers), std::move(switched_shunts)};

  return nw;
}

Network ParseNetwork(const std::string &filename) {
  std::ifstream is(filename);
  if (!is) {
    throw std::runtime_error("Failed to open file: \'" + filename + "\'");
  }
  auto nw = ParseNetwork(is);
  nw.file_name = filename;
  return nw;
}

PetscErrorCode ConvertToPS(PS ps, const Network &nw) {
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

  auto deg2rad = [](auto &&deg) { return deg * M_PI / 180.0; };

  // buses
  ps->Nbus = ps->nbus = nw.buses.size();
  ierr = PetscCalloc1(ps->Nbus, &ps->bus);
  CHKERRQ(ierr);
  std::size_t maxbusi = 0;
  for (int i = 0; i < ps->Nbus; ++i) {
    auto &dbus = ps->bus[i];
    const auto &sbus = nw.buses[i];
    dbus.bus_i = sbus.i;
    maxbusi = std::max(maxbusi, sbus.i);
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

  ps->maxbusnum = maxbusi;
  ierr = PetscCalloc1(ps->maxbusnum + 1, &ps->busext2intmap);
  CHKERRQ(ierr);
  for (int i = 0; i < ps->maxbusnum + 1; i++) {
    ps->busext2intmap[i] = -1;
  }
  for (const auto &[ext_i, int_i] : nw.bus_mapping.GetIdToIdMap()) {
    ps->busext2intmap[ext_i] = int_i;
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

    auto bus_ii = nw.bus_mapping.GetInternalIndex(dload.bus_i);
    dload.internal_i = bus_ii;
    auto &bus = ps->bus[bus_ii];
    bus.lidx[bus.nload] = i;
    bus.nload++;
  }

  // fixed_bus_shunts
  for (auto &shunt : nw.fixed_bus_shunts) {
    if (shunt.status == 0) {
      continue;
    }
    auto bus_ii = nw.bus_mapping.GetInternalIndex(shunt.i);
    if (ps->bus[bus_ii].nshunt > 0) {
      throw std::runtime_error(
          "Bus " + std::to_string(shunt.i) +
          ": more than one fixed shunt at bus not supported");
    }
    ps->bus[bus_ii].nshunt++;
    ps->bus[bus_ii].gl = shunt.gl / ps->MVAbase;
    ps->bus[bus_ii].bl = shunt.bl / ps->MVAbase;
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
    auto bus_ii = nw.bus_mapping.GetInternalIndex(dgen.bus_i);
    dgen.internal_i = bus_ii;
    auto &bus = ps->bus[bus_ii];
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
    line.internal_i = nw.bus_mapping.GetInternalIndex(line.fbus);
    line.internal_j = nw.bus_mapping.GetInternalIndex(line.tbus);
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
