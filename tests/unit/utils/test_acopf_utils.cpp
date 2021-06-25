#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libgen.h>
#include <limits>
#include <petscviewer.h>
#include <sys/stat.h>
#include <test_acopf_utils.h>

void saveToFile(Vec vector, std::string file) {

  PetscScalar *array;
  PetscInt size;

  VecGetArray(vector, &array);
  VecGetSize(vector, &size);

  std::ofstream f;

  f.open(file.c_str());
  if (f.is_open()) {
    f << size << std::endl;
    f << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < size; i++) {
      f << array[i] << std::endl;
    }
    f.close();
  } else {
    std::cout << "Failed to open file for writing Petsc Vec.\n";
  }
}

void saveToFile(double data, std::string file) {
  std::ofstream f;

  f.open(file.c_str());
  if (f.is_open()) {
    f << std::setprecision(std::numeric_limits<double>::max_digits10) << data;
    f.close();
  } else {
    std::cout << "Failed to open file for writing double.\n";
  }
}

void saveToFile(Mat a, std::string file) {
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file.c_str(), FILE_MODE_WRITE,
                        &viewer);
  MatView(a, viewer);
  PetscViewerDestroy(&viewer);
}

void readFromFile(Vec *vector, std::string file) {

  PetscErrorCode ierr;

  ierr = VecCreate(PETSC_COMM_WORLD, vector);
  CHKERRV(ierr);

  std::ifstream f(file);
  if (f.is_open()) {
    PetscInt n;

    f >> n;

    ierr = VecSetSizes(*vector, n, PETSC_DECIDE);
    CHKERRV(ierr);
    ierr = VecSetFromOptions(*vector);
    CHKERRV(ierr);

    PetscInt *indices = new PetscInt[n];
    PetscScalar *values = new PetscScalar[n];
    PetscScalar tmp;

    int i = 0;
    while (f >> tmp) {
      values[i] = tmp;
      indices[i] = i;
      i++;
    }

    f.close();

    ierr = VecSetValues(*vector, n, indices, values, INSERT_VALUES);
    CHKERRV(ierr);

    delete[] indices;
    delete[] values;
  } else {
    std::cout << "Failed to open file for reading Petsc Vec.\n";
  }
}

void readFromFile(double *data, std::string file) {
  std::ifstream f(file);
  if (f.is_open()) {
    f >> *data;
    f.close();
  } else {
    std::cout << "Failed to open file for reading double.\n";
  }
}

bool readFromFile(Mat *a, std::string file_name) {
  // See is the file exists...
  if (FILE *file = fopen(file_name.c_str(), "r")) {
    fclose(file);
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, file_name.c_str(), FILE_MODE_READ,
                          &viewer);
    MatCreate(PETSC_COMM_WORLD, a);
    MatSetType(*a, MATSEQAIJ);
    MatLoad(*a, viewer);
    PetscViewerDestroy(&viewer);
    return true;
  } else {
    return false;
  }
}

std::string getFileName(std::string file_path) {
  char *tmp = new char[file_path.size() + 1];
  std::strcpy(tmp, file_path.c_str());
  std::string res(basename(tmp));
  delete[] tmp;
  return res;
}

void validate_directory(std::string pth) {
  // std::cout << "Validating dir" << std::endl;
  DIR *dp;
  dp = opendir(pth.c_str());
  if (dp == NULL) {
    std::cout << "No dir found, creating one in place.\n";
    int status =
        mkdir(pth.c_str(), S_IRWXU | S_IRWXG | S_IROTH |
                               S_IXOTH); // Permissions to read/write/search for
                                         // group and read/search for others
    if (status != 0) {
      std::cout << "Failed to create a directory.\n";
    }
  }
  // std::cout << "Existing dir found" << std::endl;
  return;
}

