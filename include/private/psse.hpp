#pragma once

#include <array>
#include <istream>
#include <limits>
#include <string>
#include <vector>
#include <unordered_map>

#include <ps.h>

namespace exago {
namespace psse {

template <typename T>
inline constexpr T Invalid = std::numeric_limits<T>::max();

template <typename T> inline constexpr bool IsInvalid(const T &val) noexcept {
  return val == Invalid<T>;
}

struct CaseID {
  int ic{0};
  double sbase{100.0};
  int rev{34};
  int xfrrat;
  int nxfrat;
  double basfrq;
  std::array<std::string, 3> extra;
};

struct Bus {
  std::size_t i;
  std::string name{};
  double baskv{0.0};
  int ide{1};
  std::size_t area{1};
  std::size_t zone{1};
  std::size_t owner{1};
  double vm{1.0};
  double va{0.0};
  double nvhi{1.1};
  double nvlo{0.9};
  double evhi{1.1};
  double evlo{0.9};
};

struct BusRef {
  std::size_t id{Invalid<std::size_t>};
  std::string name{};
  const Bus *bus{nullptr};

  operator std::size_t() const { return id; }
  const Bus *operator->() const noexcept { return bus; }
};

struct Load {
  BusRef i;
  std::string id{"1"};
  int status{1};
  std::size_t area{Invalid<std::size_t>}; // Default is area of bus i
  std::size_t zone{Invalid<std::size_t>}; // Default is zone of bus i
  double pl{0.0};
  double ql{0.0};
  double ip{0.0};
  double iq{0.0};
  double yp{0.0};
  double yq{0.0};
  std::size_t owner{Invalid<std::size_t>}; // Default is owner of bus i
  int scale{1};
  int intrpt{0};
  double dgenp{0.0};
  double dgenq{0.0};
  int dgenm{0};
};

struct FixedBusShunt {
  BusRef i;
  std::string id{"1"};
  int status{1};
  double gl{0.0};
  double bl{0.0};
};

struct Ownership {
  std::size_t owner{0};
  double fraction{1.0};
};

struct Generator {
  BusRef i;
  std::string id{"1"};
  double pg{0.0};
  double qg{0.0};
  double qt{9999.0};
  double qb{-9999.0};
  double vs{1.0};
  BusRef ireg{0, ""};
  double mbase{Invalid<double>}; // Default is system MVA base
  double zr{0.0};
  double zx{1.0};
  double rt{0.0};
  double xt{0.0};
  double gtap{1.0};
  int stat{1};
  double rmpct{100.0};
  double pt{9999.0};
  double pb{-9999.0};
  std::array<Ownership, 4> owners{};
  int wmod{0};
  double wpf{1.0};
  std::size_t nreg{0};
};

struct Branch {
  BusRef i;
  BusRef j;
  std::string ckt;
  double r;
  double x;
  double b{0.0};
  std::string name{};
  std::array<double, 12> rates{};
  double gi{0.0};
  double bi{0.0};
  double gj{0.0};
  double bj{0.0};
  int st{1};
  int met{1};
  double len{0.0};
  std::array<Ownership, 4> owners{};
};

struct SystemSwitchingDevice {};

struct Impedence {
  double r;
  double x;
  double sbase{Invalid<double>}; // Default is system base MVA
};

struct Winding {
  /**
   * Default depends on Transformer::cw:
   * - cw in {1, 3} => 1.0
   * - cw == 2 => base voltage of corresponding bus (baskv of bus i, j or k)
   */
  double windv{Invalid<double>};
  double nomv{0.0};
  double ang{0.0};
  std::array<double, 12> rates{};
  int cod{0};
  int cont{0};
  double rma{1.1};
  double rmi{0.9};
  double vma{1.1};
  double vmi{0.9};
  int ntp{33};
  int tab{0};
  double cr{0.0};
  double cx{0.0};
  double cnxa{0.0};
  std::size_t node{0};
};

struct Transformer {
  BusRef i;
  BusRef j;
  BusRef k;
  std::string ckt{"1"};
  int cw{1};
  int cz{1};
  int cm{1};
  double mag1{0.0};
  double mag2{0.0};
  int nmetr{2};
  std::string name{};
  int stat{1};
  std::array<Ownership, 4> owners{};
  std::string vecgrp{};
  Impedence imp12;
  Impedence imp23;
  Impedence imp31;
  double vmstar{1.0};
  double anstar{0.0};
  std::array<Winding, 3> windings{};
};

struct AreaInterchange {
  std::size_t i;
  BusRef isw;
  double pdes{0.0};
  double ptol{10.0};
  std::string arname{};
};

struct TwoTerminalDCLine {};

struct VSCDCLine {};

struct ImpedanceCorrection {};

struct MultiTerminalDCLine {};

struct MultiSectionLineGroup {};

struct Zone {
  std::size_t i;
  std::string zoname{};
};

struct InterAreaTransfer {};

struct Owner {
  std::size_t i;
  std::string owname{};
};

struct FACTSDevice {};

struct SwitchedShunt {
  struct Block {
    std::size_t n{0};
    double b{0.0};
  };

  BusRef i;
  int modsw{1};
  int adjm{0};
  int stat{1};
  double vswhi{1.0};
  double vswlo{1.0};
  BusRef swreg{0, ""};
  double rmpct{100.0};
  std::string rmidnt{};
  double binit{0.0};
  std::array<Block, 8> blocks{};
  BusRef nreg{0, ""};
};

struct GNEDevice {};

struct InductionMachine {};

class BusMapping {
public:
  struct Optional {
    operator bool() { return value; }
    bool value{false};
  };

  BusMapping(std::vector<Bus> &buses, bool require_unique_names);

  bool RequireUniqueNames() const noexcept { return req_unique_names_; }

  bool HasBus(std::size_t bus_number) const;
  bool HasBus(const std::string &bus_name) const;

  std::size_t GetInternalIndex(std::size_t bus_number) const;
  std::size_t GetInternalIndex(const std::string &bus_name) const;

  std::size_t GetBusNumber(const std::string &bus_name) const;
  const std::string &GetBusName(std::size_t bus_number) const;

  const Bus &GetBus(const std::string &bus_name) const;
  const Bus &GetBus(std::size_t bus_number) const;

  void Resolve(BusRef &bus, Optional opt = Optional{false}) const;

  const auto &GetIdToIdMap() const { return id_map_; }
  const auto &GetNameToIdMap() const { return name_map_; }

private:
  bool req_unique_names_;
  std::vector<Bus> &buses_;
  std::unordered_map<std::size_t, std::size_t> id_map_;
  std::unordered_map<std::string, std::size_t> name_map_;
};

struct Network {
  Network(bool ref_bus_names, CaseID &&, std::vector<Bus> &&,
          std::vector<Load> &&, std::vector<FixedBusShunt> &&,
          std::vector<Generator> &&, std::vector<Branch> &&,
          std::vector<Transformer> &&, std::vector<AreaInterchange> &&,
          std::vector<Zone> &&, std::vector<Owner> &&,
          std::vector<SwitchedShunt> &&);

  void ResolveBusIds();
  void ResolveDefaults();

  std::string file_name;
  CaseID case_id;
  std::vector<Bus> buses;
  BusMapping bus_mapping;
  std::vector<Load> loads;
  std::vector<FixedBusShunt> fixed_bus_shunts;
  std::vector<Generator> generators;
  std::vector<Branch> branches;
  // std::vector<SystemSwitchingDevice> system_switching_devices;
  std::vector<Transformer> transformers;
  std::vector<AreaInterchange> area_interchanges;
  // std::vector<TwoTerminalDCLine> two_terminal_dc_lines;
  // std::vector<VSCDCLine> vsc_dc_lines;
  // std::vector<ImpedanceCorrection> impedance_corrections;
  // std::vector<MultiTerminalDCLine> multi_terminal_dc_lines;
  // std::vector<MultiSectionLineGroup> multi_section_line_groups;
  std::vector<Zone> zones;
  // std::vector<InterAreaTransfer> inter_area_transfers;
  std::vector<Owner> owners;
  // std::vector<FACTSDevice> facts_devices;
  std::vector<SwitchedShunt> switched_shunts;
  // std::vector<GNEDevice> gne_devices;
};

Network ParseNetwork(std::istream &is);
Network ParseNetwork(const std::string &filename);

PetscErrorCode ConvertToPS(PS ps, const Network &nw);

} // namespace psse
} // namespace exago
