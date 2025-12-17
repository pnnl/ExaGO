#pragma once

#include <array>
#include <istream>
#include <string>
#include <vector>

#include <ps.h>

namespace exago {
namespace psse {

struct CaseID {
  int ic;
  double sbase;
  int rev;
  int xfrrat;
  int nxfrat;
  double basfrq;
  std::array<std::string, 3> extra;
};

struct Bus {
  std::size_t i;
  std::string name;
  double baskv;
  int ide;
  int area;
  int zone;
  int owner;
  double vm;
  double va;
  double nvhi{1.1};
  double nvlo{0.9};
  double evhi{1.1};
  double evlo{0.9};
};

struct Load {
  std::size_t i;
  std::string i_bus_name;
  std::string id;
  int status;
  int area;
  int zone;
  double pl;
  double ql;
  double ip;
  double iq;
  double yp;
  double yq;
  int owner;
  int scale{1};
  int intrpt{0};
};

struct FixedBusShunt {
  std::size_t i;
  std::string i_bus_name;
  std::string id;
  int status;
  double gl;
  double bl;
};

struct Ownership {
  int owner;
  double fraction;
};

struct Generator {
  std::size_t i;
  std::string i_bus_name;
  std::string id;
  double pg;
  double qg;
  double qt;
  double qb;
  double vs;
  std::size_t ireg;
  std::string ireg_bus_name;
  double mbase;
  double zr;
  double zx;
  double rt;
  double xt;
  double gtap;
  int stat;
  double rmpct;
  double pt;
  double pb;
  std::array<Ownership, 4> owners;
  int wmod{0};
  double wpf{1.0};
};

struct Branch {
  std::size_t i;
  std::string i_bus_name;
  std::size_t j;
  std::string j_bus_name;
  std::string ckt;
  double r;
  double x;
  double b;
  double ratea;
  double rateb;
  double ratec;
  double gi;
  double bi;
  double gj;
  double bj;
  int st;
  int met;
  double len;
  std::array<Ownership, 4> owners;
};

struct Impedence {
  double r;
  double x;
  double sbase;
};

struct Winding {
  double windv;
  double nomv;
  double ang;
  double rata;
  double ratb;
  double ratc;
  int cod;
  int cont;
  double rma;
  double rmi;
  double vma;
  double vmi;
  int ntp;
  int tab;
  double cr;
  double cx;
  double cnxa;
};

struct Transformer {
  std::size_t i;
  std::string i_bus_name;
  std::size_t j;
  std::string j_bus_name;
  std::size_t k;
  std::string k_bus_name;
  std::string ckt;
  int cw;
  int cz;
  int cm;
  double mag1;
  double mag2;
  int nmetr;
  std::string name;
  int stat;
  std::array<Ownership, 4> owners;
  Impedence imp12;
  Impedence imp23;
  Impedence imp31;
  double vmstar;
  double anstar;
  std::array<Winding, 3> windings;
};

struct Area {
  int i;
  std::size_t isw;
  std::string isw_bus_name; // TODO: resolve
  double pdes;
  double ptol;
  std::string arname;
};

class BusMapping {
public:
  struct Optional {
    operator bool() { return value; }
    bool value{false};
  };

  BusMapping(const std::vector<Bus> &buses);

  std::size_t getInternalIndex(std::size_t bus_number) const;
  std::size_t getInternalIndex(const std::string &bus_name) const;

  std::size_t getBusNumber(const std::string &bus_name) const;
  const std::string &getBusName(std::size_t bus_number) const;

  const Bus &getBus(const std::string &bus_name) const;
  const Bus &getBus(std::size_t bus_number) const;

  void resolve(std::size_t &bus_number, std::string &bus_name,
               Optional opt = Optional{false}) const;

private:
  const std::vector<Bus> &buses_;
  std::unordered_map<std::size_t, std::size_t> id_map_;
  std::unordered_map<std::string, std::size_t> name_map_;
};

struct Network {
  Network(CaseID &&, std::vector<Bus> &&, std::vector<Load> &&,
          std::vector<FixedBusShunt> &&, std::vector<Generator> &&,
          std::vector<Branch> &&, std::vector<Transformer> &&);

  void resolveBusIds();

  std::string file_name;
  CaseID case_id;
  std::vector<Bus> buses;
  BusMapping bus_mapping;
  std::vector<Load> loads;
  std::vector<FixedBusShunt> shunts;
  std::vector<Generator> generators;
  std::vector<Branch> branches;
  std::vector<Transformer> transformers;
};

Network parse_network(std::istream &is);
Network parse_network(const std::string &filename);

PetscErrorCode convert_to_ps(PS ps, const Network &nw);

} // namespace psse
} // namespace exago
