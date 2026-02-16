#include <algorithm>
#include <cctype>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <variant>

#include <psse.hpp>
#include <psimpl.h>
#include <utils.h>

namespace exago {
namespace psse {

template <typename T> bool Approx(T a, T b, T rtol = 1.0e-5) {
  return std::fabs(a - b) <= rtol * std::max(std::fabs(a), std::fabs(b));
}

[[noreturn]] void Error(const std::string &message) {
  ExaGOLog(EXAGO_LOG_ERROR, message);
  throw ExaGOError(message);
}

void Warn(const std::string &message) { ExaGOLog(EXAGO_LOG_WARN, message); }

std::string Strip(std::string str) {
  auto notspace = [](char c) { return !std::isspace(c); };
  str.erase(begin(str), std::find_if(begin(str), end(str), notspace));
  str.erase(std::find_if(rbegin(str), rend(str), notspace).base(), end(str));
  return str;
}

int CountWords(const std::string &str) {
  std::stringstream stream(str);
  std::istream_iterator<std::string> it(stream);
  std::istream_iterator<std::string> end;
  return static_cast<int>(std::distance(it, end));
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
      Error("First character expected to be single or double quote");
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
  bool IsTypeInt() const { return std::holds_alternative<int>(var_); }

  int GetInt() const {
    if (IsTypeInt()) {
      return std::get<int>(var_);
    } else {
      return 0;
    }
  }

  bool IsTypeString() const {
    return std::holds_alternative<std::string>(var_);
  }

  std::string GetString() const {
    if (IsTypeString()) {
      return std::get<std::string>(var_);
    } else {
      return "";
    }
  }

  BusRef ToBusRef() const {
    return BusRef{static_cast<std::size_t>(GetInt()), GetString()};
  }

private:
  friend std::istream &operator>>(std::istream &is, IntOrStringParse &p) {
    SkipLeadingWhitespace(is);

    // Check for opening quote
    auto qc = is.peek();
    if (qc == '\'' || qc == '\"') {
      QuoteStringParse qs;
      is >> qs;
      p.var_ = Strip(qs);
    } else {
      is >> p.var_.emplace<int>();
    }

    return is;
  }

  std::variant<int, std::string> var_;
};

class LineItemStream;

class CheckedLineItemStream {
public:
  CheckedLineItemStream(LineItemStream &lis) : lis_(lis) {}

  template <typename T> CheckedLineItemStream &operator>>(T &item);

private:
  LineItemStream &lis_;
};

struct CheckStream {};

inline constexpr CheckStream Check{};

class LineItemStream : public std::istream {
public:
  LineItemStream() = delete;
  LineItemStream(std::istream &is) : std::istream(is.rdbuf()), is_(&is) {
    do {
      NextLine();
    } while (StartsWith("@!"));
  }

  CheckedLineItemStream Checked() { return CheckedLineItemStream(*this); }

  operator bool() const { return !q_.empty(); }

  std::size_t Size() const noexcept { return q_.size(); }
  const std::string &Raw() const noexcept { return line_; }

  bool StartsWith(const std::string &sub) const {
    return q_.front().find(sub) == 0;
  }

  LineItemStream &NextLine() {
    line_ = ReadLine(*is_);
    std::istringstream iss(line_);
    q_.clear();
    for (std::string item; std::getline(iss, item, ',');) {
      q_.push_back(Strip(item));
    }
    return *this;
  }

  std::string String() const { return line_; }

  template <typename T> LineItemStream &operator>>(T &item) {
    if (!q_.front().empty()) {
      std::istringstream(q_.front()) >> item;
    }
    q_.pop_front();
    return *this;
  }

  CheckedLineItemStream operator>>(CheckStream) { return Checked(); }

private:
  friend std::string ReadLine(LineItemStream &in) {
    std::string out;
    for (auto &&item : in.q_) {
      out += " " + item;
    }
    return out;
  }

  std::istream *is_{nullptr};

  std::string line_;
  std::deque<std::string> q_;
};

template <typename T>
inline CheckedLineItemStream &CheckedLineItemStream::operator>>(T &item) {
  if (lis_.Size() > 0) {
    lis_ >> item;
  }
  return *this;
}

CaseID ParseCaseID(std::istream &is) {
  CaseID cid;
  std::string record;
  std::getline(is, record, '/');
  std::istringstream rec_stream(record);
  LineItemStream lis(rec_stream);
  lis >> cid.ic >> cid.sbase >> cid.rev;
  int supported[] = {32, 33, 34};
  if (std::find(std::begin(supported), std::end(supported), cid.rev) ==
      std::end(supported)) {
    Error("PSS(R)E Power Flow Raw Data Files are supported only for "
          "versions 32, 33 and 34. This file is tagged version " +
          std::to_string(cid.rev));
  }
  lis >> cid.xfrrat >> cid.nxfrat >> cid.basfrq;

  cid.extra[0] = ReadLine(is);
  cid.extra[1] = ReadLine(is);
  cid.extra[2] = ReadLine(is);
  return cid;
}

void SkipSystemWideData(std::istream &is) {
  LineItemStream lis(is);
  while (!lis.StartsWith("0 /")) {
    lis.NextLine();
  }
}

struct Parser {
  int rev{34};
  bool ref_bus_names{false};

  void ParseRecord(LineItemStream &lis, Bus &bus) {
    QuoteStringParse name;
    lis >> bus.i;
    lis &&lis >> name;
    bus.name = Strip(name);
    if (bus.name.empty()) {
      bus.name = "BUS" + std::to_string(bus.i);
    }
    auto clis = lis.Checked();
    clis >> bus.baskv >> bus.ide >> bus.area >> bus.zone >> bus.owner >>
        bus.vm >> bus.va >> bus.nvhi >> bus.nvlo >> bus.evhi >> bus.evlo;
  }

  void ParseRecord(LineItemStream &lis, Load &ld) {
    IntOrStringParse i;
    QuoteStringParse id;
    lis >> i;
    if (i.IsTypeString()) {
      ref_bus_names = true;
    }
    ld.i = i.ToBusRef();
    auto clis = lis.Checked();
    clis >> id;
    ld.id = Strip(id);
    clis >> ld.status >> ld.area >> ld.zone >> ld.pl >> ld.ql >> ld.ip >>
        ld.iq >> ld.yp >> ld.yq >> ld.owner >> ld.scale >> ld.intrpt;
    if (rev >= 34) {
      clis >> ld.dgenp >> ld.dgenq >> ld.dgenm;
    }
  }

  void ParseRecord(LineItemStream &lis, FixedBusShunt &sh) {
    IntOrStringParse i;
    QuoteStringParse id;
    lis >> i;
    if (i.IsTypeString()) {
      ref_bus_names = true;
    }
    sh.i = i.ToBusRef();
    auto clis = lis.Checked();
    clis >> id;
    sh.id = Strip(id);
    clis >> sh.status >> sh.gl >> sh.bl;
  }

  void ParseRecord(LineItemStream &lis, Generator &gen) {
    IntOrStringParse i;
    IntOrStringParse ireg;
    QuoteStringParse id;
    lis >> i;
    gen.i = i.ToBusRef();
    auto clis = lis.Checked();
    clis >> id;
    gen.id = Strip(id);
    clis >> gen.pg >> gen.qg >> gen.qt >> gen.qb >> gen.vs >> ireg;
    if (i.IsTypeString() || ireg.IsTypeString()) {
      ref_bus_names = true;
    }
    gen.ireg = ireg.ToBusRef();
    clis >> gen.mbase >> gen.zr >> gen.zx >> gen.rt >> gen.xt >> gen.gtap >>
        gen.stat >> gen.rmpct >> gen.pt >> gen.pb;
    for (std::size_t oi = 0; oi < 4; ++oi) {
      clis >> gen.owners[oi].owner >> gen.owners[oi].fraction;
    }
    clis >> gen.wmod >> gen.wpf;
    if (rev >= 34) {
      clis >> gen.nreg;
    }
  }

  void ParseRecord(LineItemStream &lis, Branch &br) {
    IntOrStringParse i;
    IntOrStringParse j;
    QuoteStringParse ckt;
    QuoteStringParse name;
    lis >> i >> j >> ckt >> br.r >> br.x;
    if (i.IsTypeString() || j.IsTypeString()) {
      ref_bus_names = true;
    }
    br.i = i.ToBusRef();
    br.j = j.ToBusRef();
    br.ckt = Strip(ckt);
    auto clis = lis.Checked();
    clis >> br.b;
    if (rev >= 34) {
      clis >> name;
      br.name = Strip(name);
    }
    clis >> br.rates[0] >> br.rates[1] >> br.rates[2];
    if (rev >= 34) {
      for (std::size_t ri = 3; ri < br.rates.size(); ++ri) {
        clis >> br.rates[ri];
      }
    }
    clis >> br.gi >> br.bi >> br.gj >> br.bj >> br.st >> br.met >> br.len;
    for (std::size_t oi = 0; oi < 4; ++oi) {
      clis >> br.owners[oi].owner >> br.owners[oi].fraction;
    }
  }

  void ParseRecord(LineItemStream &, SystemSwitchingDevice &) {}

  Winding ParseTransformerWinding(LineItemStream &lis) {
    Winding w;
    auto clis = lis.Checked();
    clis >> w.windv >> w.nomv >> w.ang >> w.rates[0] >> w.rates[1] >>
        w.rates[2];
    if (rev >= 34) {
      clis >> w.rates[3] >> w.rates[4] >> w.rates[5] >> w.rates[6] >>
          w.rates[7] >> w.rates[8] >> w.rates[9] >> w.rates[10] >> w.rates[11];
    }
    clis >> w.cod >> w.cont >> w.rma >> w.rmi >> w.vma >> w.vmi >> w.ntp >>
        w.tab >> w.cr >> w.cx >> w.cnxa;
    return w;
  }

  void ParseRecord(LineItemStream &lis, Transformer &tr) {
    IntOrStringParse i;
    IntOrStringParse j;
    IntOrStringParse k;
    QuoteStringParse ckt;
    QuoteStringParse name;
    lis >> i >> j >> k;
    if (i.IsTypeString() || j.IsTypeString() || k.IsTypeString()) {
      ref_bus_names = true;
    }
    tr.i = i.ToBusRef();
    tr.j = j.ToBusRef();
    tr.k = k.ToBusRef();
    auto clis = lis.Checked();
    clis >> ckt;
    tr.ckt = Strip(ckt);
    clis >> tr.cw >> tr.cz >> tr.cm >> tr.mag1 >> tr.mag2 >> tr.nmetr;
    clis >> name;
    tr.name = Strip(name);
    clis >> tr.stat;
    for (std::size_t oi = 0; oi < 4; ++oi) {
      clis >> tr.owners[0].owner >> tr.owners[0].fraction;
    }
    if (tr.k == 0) {
      // two-winding (3 more rows)
      lis.NextLine() >> tr.imp12.r >> tr.imp12.x >> Check >> tr.imp12.sbase;
      tr.windings[0] = ParseTransformerWinding(lis.NextLine());
      lis.NextLine() >> Check >> tr.windings[1].windv >> tr.windings[1].nomv;
    } else {
      // three-winding (4 more rows)
      lis.NextLine() >> tr.imp12.r >> tr.imp12.x >> Check >> tr.imp12.sbase >>
          tr.imp23.r >> tr.imp23.x >> tr.imp23.sbase >> tr.imp31.r >>
          tr.imp31.x >> tr.imp31.sbase >> Check >> tr.vmstar >> tr.anstar;
      tr.windings[0] = ParseTransformerWinding(lis.NextLine());
      tr.windings[1] = ParseTransformerWinding(lis.NextLine());
      tr.windings[2] = ParseTransformerWinding(lis.NextLine());
    }
  }

  void ParseRecord(LineItemStream &lis, AreaInterchange &area) {
    IntOrStringParse isw;
    QuoteStringParse arname;
    lis >> area.i;
    auto clis = lis.Checked();
    clis >> isw;
    if (isw.IsTypeString()) {
      ref_bus_names = true;
    }
    area.isw = isw.ToBusRef();
    clis >> area.pdes >> area.ptol >> arname;
    area.arname = Strip(arname);
  }

  void ParseRecord(LineItemStream &, TwoTerminalDCLine &) {}

  void ParseRecord(LineItemStream &, VSCDCLine &) {}

  void ParseRecord(LineItemStream &, ImpedanceCorrection &) {}

  void ParseRecord(LineItemStream &, MultiTerminalDCLine &) {}

  void ParseRecord(LineItemStream &, MultiSectionLineGroup &) {}

  void ParseRecord(LineItemStream &lis, Zone &zone) {
    QuoteStringParse zoname;
    lis >> zone.i >> Check >> zoname;
    zone.zoname = Strip(zoname);
  }

  void ParseRecord(LineItemStream &, InterAreaTransfer &) {}

  void ParseRecord(LineItemStream &lis, Owner &owner) {
    QuoteStringParse owname;
    lis >> owner.i >> Check >> owname;
    owner.owname = Strip(owname);
  }

  void ParseRecord(LineItemStream &, FACTSDevice &) {}

  void ParseRecord(LineItemStream &lis, SwitchedShunt &sh) {
    IntOrStringParse i;
    IntOrStringParse swreg;
    QuoteStringParse rmidnt;
    lis >> i;
    sh.i = i.ToBusRef();
    auto clis = lis.Checked();
    clis >> sh.modsw >> sh.adjm >> sh.stat >> sh.vswhi >> sh.vswlo >> swreg;
    sh.swreg = swreg.ToBusRef();
    if (i.IsTypeString() || swreg.IsTypeString()) {
      ref_bus_names = true;
    }
    clis >> sh.rmpct >> rmidnt;
    sh.rmidnt = Strip(rmidnt);
    clis >> sh.binit;
    for (std::size_t bi = 0; bi < 8; ++bi) {
      clis >> sh.blocks[bi].n >> sh.blocks[bi].b;
    }
    if (rev >= 34 && lis) {
      IntOrStringParse nreg;
      clis >> nreg;
      if (nreg.IsTypeString()) {
        ref_bus_names = true;
      }
      sh.nreg = nreg.ToBusRef();
    }
  }

  void ParseRecord(LineItemStream &, GNEDevice &) {}
  void ParseRecord(LineItemStream &, InductionMachine &) {}

  template <typename T> std::vector<T> ParseRecords(std::istream &is) {
    std::vector<T> recs;
    while (is && is.peek() != 'Q') {
      LineItemStream lis(is);
      if (lis.StartsWith("0 /")) {
        break;
      }
      auto &rec = recs.emplace_back();
      ParseRecord(lis, rec);
    }
    return recs;
  }
};

BusMapping::BusMapping(std::vector<Bus> &buses, bool require_unique_names)
    : req_unique_names_(require_unique_names), buses_(buses) {
  for (std::size_t i = 0; i < buses_.size(); ++i) {
    auto [_, ins] = id_map_.emplace(buses_[i].i, i);
    if (!ins) {
      Error("PSS(R)E parser: Encountered duplicate bus number: " +
            std::to_string(buses_[i].i));
    }
  }
  if (req_unique_names_) {
    for (std::size_t i = 0; i < buses_.size(); ++i) {
      auto [_, ins] = name_map_.emplace(buses_[i].name, i);
      if (!ins) {
        Error("PSS(R)E parser: Encountered duplicate bus name: " +
              buses_[i].name);
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

void BusMapping::Resolve(BusRef &busref, BusMapping::Optional optional) const {
  if (req_unique_names_) {
    if (busref.id == 0 && busref.name.empty()) {
      if (optional) {
        return;
      } else {
        Error("Cannot resolve bus id");
      }
    }
    if (busref.id == 0) {
      if (!HasBus(busref.name)) {
        Error("PSS(R)E parser: Bus \'" + busref.name + "\' does not exist");
      }
      const Bus &bus = GetBus(busref.name);
      busref.id = bus.i;
      busref.bus = &bus;
    }
    if (busref.name.empty()) {
      if (!HasBus(busref.id)) {
        Error("PSS(R)E parser: Bus " + std::to_string(busref.id) +
              " does not exist");
      }
      const Bus &bus = GetBus(busref.id);
      busref.name = bus.name;
      busref.bus = &bus;
    }
  } else {
    if (busref.id == 0) {
      if (optional) {
        return;
      }
      Error("PSS(R)E parser: Bus reference number must be positive, got 0");
    }
    busref.bus = &GetBus(busref.id);
  }
}

void Network::ResolveBusIds() {
  using Optional = BusMapping::Optional;
  for (auto &load : loads) {
    bus_mapping.Resolve(load.i);
  }
  for (auto &shunt : fixed_bus_shunts) {
    bus_mapping.Resolve(shunt.i);
  }
  for (auto &gen : generators) {
    bus_mapping.Resolve(gen.i);
    bus_mapping.Resolve(gen.ireg, Optional{true});
  }
  for (auto &br : branches) {
    bus_mapping.Resolve(br.i);
    bus_mapping.Resolve(br.j);
  }
  for (auto &tr : transformers) {
    bus_mapping.Resolve(tr.i);
    bus_mapping.Resolve(tr.j);
    bus_mapping.Resolve(tr.k, Optional{true});
  }
  for (auto &ar : area_interchanges) {
    bus_mapping.Resolve(ar.isw, Optional{true});
  }
  for (auto &sh : switched_shunts) {
    bus_mapping.Resolve(sh.i);
    bus_mapping.Resolve(sh.swreg, Optional{true});
    bus_mapping.Resolve(sh.nreg, Optional{true});
  }
}

void Network::ResolveDefaults() {
  for (auto &load : loads) {
    if (load.i) {
      if (IsInvalid(load.area)) {
        load.area = load.i->area;
      }
      if (IsInvalid(load.zone)) {
        load.zone = load.i->zone;
      }
      if (IsInvalid(load.owner)) {
        load.owner = load.i->owner;
      }
    }
  }
  auto resolve_o1 = [](auto &comp) {
    if (comp.owners[0].owner == 0) {
      if (comp.i) {
        comp.owners[0].owner = comp.i->owner;
      }
    }
  };
  for (auto &gen : generators) {
    if (IsInvalid(gen.mbase)) {
      gen.mbase = case_id.sbase;
    }
    resolve_o1(gen);
  }
  for (auto &br : branches) {
    resolve_o1(br);
  }
  for (auto &tr : transformers) {
    resolve_o1(tr);
    for (std::size_t wi = 0; wi < 3; ++wi) {
      auto &w = tr.windings[wi];
      if (IsInvalid(w.windv)) {
        if (tr.cw == 2) {
          auto &busref = wi == 0 ? tr.i : wi == 1 ? tr.j : tr.k;
          if (busref) {
            w.windv = busref->baskv;
          }
        } else {
          w.windv = 1.0;
        }
      }
    }
  }
}

Network::Network(bool ref_bus_names, CaseID &&cid, std::vector<Bus> &&bus,
                 std::vector<Load> &&load, std::vector<FixedBusShunt> &&fbshunt,
                 std::vector<Generator> &&gen, std::vector<Branch> &&branch,
                 std::vector<Transformer> &&trans,
                 std::vector<AreaInterchange> &&area, std::vector<Zone> &&zone,
                 std::vector<Owner> &&owner,
                 std::vector<SwitchedShunt> &&swshunt)
    : case_id(std::move(cid)), buses(std::move(bus)),
      bus_mapping(buses, ref_bus_names), loads(std::move(load)),
      fixed_bus_shunts(std::move(fbshunt)), generators(std::move(gen)),
      branches(std::move(branch)), transformers(std::move(trans)),
      area_interchanges(std::move(area)), zones(std::move(zone)),
      owners(std::move(owner)), switched_shunts(std::move(swshunt)) {
  ResolveBusIds();
  ResolveDefaults();
}

Network ParseNetwork(std::istream &is) {
  auto case_id = ParseCaseID(is);
  if (case_id.rev >= 34) {
    SkipSystemWideData(is);
  }
  auto parser = Parser{case_id.rev};
  auto buses = parser.ParseRecords<Bus>(is);
  auto loads = parser.ParseRecords<Load>(is);
  auto fixed_bus_shunts = parser.ParseRecords<FixedBusShunt>(is);
  auto generators = parser.ParseRecords<Generator>(is);
  auto branches = parser.ParseRecords<Branch>(is);
  // { TODO
  if (case_id.rev >= 34) {
    auto system_switching_devices =
        parser.ParseRecords<SystemSwitchingDevice>(is);
  }
  // }
  auto transformers = parser.ParseRecords<Transformer>(is);
  auto area_interchanges = parser.ParseRecords<AreaInterchange>(is);
  // { TODO
  auto two_terminal_dc_lines = parser.ParseRecords<TwoTerminalDCLine>(is);
  auto vsc_dc_lines = parser.ParseRecords<VSCDCLine>(is);
  auto impedance_corrections = parser.ParseRecords<ImpedanceCorrection>(is);
  auto multi_terminal_dc_lines = parser.ParseRecords<MultiTerminalDCLine>(is);
  auto multi_section_line_groups =
      parser.ParseRecords<MultiSectionLineGroup>(is);
  // }
  auto zones = parser.ParseRecords<Zone>(is);
  // { TODO
  auto inter_area_transfers = parser.ParseRecords<InterAreaTransfer>(is);
  // }
  auto owners = parser.ParseRecords<Owner>(is);
  // { TODO
  auto facts_devices = parser.ParseRecords<FACTSDevice>(is);
  // }
  auto switched_shunts = parser.ParseRecords<SwitchedShunt>(is);

  // { TODO
  auto gne_devices = parser.ParseRecords<GNEDevice>(is);
  auto induction_machines = parser.ParseRecords<InductionMachine>(is);
  // }

  Network nw{parser.ref_bus_names,
             std::move(case_id),
             std::move(buses),
             std::move(loads),
             std::move(fixed_bus_shunts),
             std::move(generators),
             std::move(branches),
             std::move(transformers),
             std::move(area_interchanges),
             std::move(zones),
             std::move(owners),
             std::move(switched_shunts)};

  return nw;
}

Network ParseNetwork(const std::string &filename) {
  std::ifstream is(filename);
  auto quoted_filename = "\'" + filename + "\'";
  if (!is) {
    Error("Failed to open file: " + quoted_filename);
  }
  ExaGOLog(EXAGO_LOG_INFO, "Parsing PSS(R)E raw data file: " + quoted_filename);
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
    dbus.Vmax = 1.1;
    dbus.Vmin = 0.9;
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
      Warn("Bus " + std::to_string(shunt.i.id) +
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
      if (!Approx(dgen.vs, bus.vm)) {
        std::stringstream ss;
        ss << "Generator at bus " << bus.bus_i << std::fixed
           << ": voltage setpoint (" << dgen.vs
           << ") different from bus voltage magnitude (" << bus.vm << ")";
        Warn(ss.str());
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
    dline.rateA = sline.rates[0];
    dline.rateB = sline.rates[1];
    dline.rateC = sline.rates[2];
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
    dline.rateA = sline.windings[0].rates[0];
    dline.rateB = sline.windings[0].rates[1];
    dline.rateC = sline.windings[0].rates[2];
    // the rest are skipped

    configure_line(dline);
  }

  auto unsupported_contents_error = [&nw](const std::string &label) {
    std::stringstream ss;
    ss << "File \'" << nw.file_name << "\' contains \'" << label
       << "\' data, which currently cannot be "
          "represented in the PS data structure";
    ExaGOLog(EXAGO_LOG_WARN, ss.str());
  };

  if (!nw.area_interchanges.empty()) {
    unsupported_contents_error("AREA");
  }
  if (!nw.zones.empty()) {
    unsupported_contents_error("ZONE");
  }
  if (!nw.owners.empty()) {
    unsupported_contents_error("OWNER");
  }
  if (!nw.switched_shunts.empty()) {
    unsupported_contents_error("SWITCHED SHUNT");
  }

  PetscFunctionReturn(0);
}

} // namespace psse
} // namespace exago
