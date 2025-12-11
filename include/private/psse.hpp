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
  int i;
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
  int i; // TODO: could also be specified as bus name (???)
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
  int i; // TODO: could be bus name (???)
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
  int i; // TODO: could be bus name (???)
  std::string id;
  double pg;
  double qg;
  double qt;
  double qb;
  double vs;
  int ireg; // TODO: could be bus name (???)
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
  int i; // TODO: could be bus name (???)
  int j; // TODO: could be bus name (???)
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
  int i;
  int j;
  int k;
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
  int isw; // TODO: could be bus name (???)
  double pdes;
  double ptol;
  std::string arname;
};

struct Network {
  std::string file_name;
  CaseID case_id;
  std::vector<Bus> buses;
  std::unordered_map<int, int> bus_id_map;
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
