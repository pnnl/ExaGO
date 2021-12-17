#include <petsc.h>
#include <opflow.h>
#include <utils.h>

/**
 * We don't have a great way to verify that the output from this is what we
 * expect from only within the test. We have to rely on some external validation
 * that the messages logged with ExaGOLog are logged on every rank, and that
 * messages from ExaGOLogSingleRank are only logged on the root rank.
 */
static char appname[] = "test_logger";
int main(int argc, char **argv) {
  MPI_Comm comm = MPI_COMM_WORLD;
  PetscErrorCode ierr;
  ExaGOInitialize(comm, nullptr, nullptr, appname, nullptr);

  // Set the minimum log level to 0 so all messages are logged
  ierr = ExaGOLogSetMinLogLevel(0);
  ExaGOCheckError(ierr);

  // The corresponding script will search for XXX and ensure this message is
  // only logged one time per run, no matter how many processes it's launched
  // with
  ExaGOLog(comm, 10, "XXX: {}",
           "This message should only be printed on a single rank");

  // The same script will ensure this message occurs as many times as processes
  // the executable is launched with
  ExaGOLog(10, "YYY: {}", "This message should be printed on every rank");

  ExaGOFinalize();
  return 0;
}
