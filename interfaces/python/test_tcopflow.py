import exago


def initialize(appname):
    exago.initialize(appname)


def finalize():
    exago.finalize()


class TCOPFLOW:
    def __init__(self):
        self.tcopflow = None

    def create(self, netfile, ploadprofile=None, qloadprofile=None, windgenprofile=None):
        self.tcopflow = exago.TCOPFLOW()
        self.tcopflow.set_network_data(netfile)

        if ploadprofile is not None or qloadprofile is not None:
            self.tcopflow.set_load_profiles(ploadprofile, qloadprofile)

        if windgenprofile is not None:
            self.tcopflow.set_wind_gen_profiles(windgenprofile)

    def solve(self):
        self.tcopflow.solve()

    def print_solution(self, index=0):
        self.tcopflow.print_solution(index)

    def save_solution(self, file_name):
        self.tcopflow.save_solution_all(exago.MATPOWER, file_name)


initialize("tcopflow")
tcopflow = TCOPFLOW()
tcopflow.create("datafiles/case9/case9mod.m")
tcopflow.solve()
tcopflow.print_solution()
tcopflow.save_solution("tcopflowout")
del tcopflow
finalize()