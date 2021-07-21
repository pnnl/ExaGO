#pragma once

#include <string>
#include <common.h>

/**
 * @brief Writes a PetscVec to _file_ in ASCII format.
 *
 * @param[in] vector PetscVec to write to _file_.
 * @param[in] file The absolute file path to be written to.
 *
 * @pre _vector_ is a constructed and populated vector.
 *
 * @post _file_ has the number of elements in the vector as the first
 *              entry in the file. The values of the vector follow.
 * @post _file_ contains data with an appropriate number of sig. figs.
 */
void saveToFile(Vec vector, std::string file);

/**
 * @brief Writes a double to _file_ in ASCII format.
 *
 * @param[in] data a double to write to _file_.
 * @param[in] file The absolute file path to be written to.
 *
 * @post _data_ is stored in file to 16 sig. fig.
 */
void saveToFile(double data, std::string file);

/**
 * @brief Writes a Petsc Matrix to _file_ in ASCII format.
 *
 * @param[in] a a Petsc Matrix to write to _file_.
 * @param[in] file The absolute file path to be written to.
 *
 * @post _a_ is stored in file in binary format.
 */
void saveToFile(Mat a, std::string file);

/**
 * @brief Reads a PetscVec from _file_ in ASCII format.
 *
 * @param[out] vector PetscVec to be read from _file_.
 * @param[in] file The absolute file path to be written to.
 *
 * @pre _file_ contains the number of elements to be stored into vec as the first
 *             entry in the file. The other values of the vector follow.
 * @pre _file_ stores the data in enough precision to be accurately stored by _vector_.
 *
 * @post _vector_ is filled with data from file.
 */

void readFromFile(Vec* vector, std::string file);

/**
 * @brief Reads a double from _file_ in ASCII format.
 *
 * @param[out] data double that is assigned data from _file_.
 * @param[in] file The absolute file path to be written to.
 *
 * @post _data_ is filled with data from file to 16 sig. figs.
 */

void readFromFile(double* data, std::string file);

/**
 * @brief Reads a Petsc Matrix from _file_ in binary format.
 *
 * @param[out] data double that is assigned data from _file_.
 * @param[in] file The absolute file path to be written to.
 * 
 * @return true if the matrix was read succesfully, false if not
 *
 * @post _a_ is filled with from a binary file.
 * @post if _file_ does not exist, then the function will silently fail.
 *       Since some answers do not store inequality constraints, this was
 *       considered desirable behaviour at the time of writing.
 */

bool readFromFile(Mat *a, std::string file);

/**
 * @brief Strips the filename from an absoulte file path.
 *
 * @param[in] file_path the absoulte file path of the file in question.
 *
 * @return The stripped filename including the file extension.
 *
 * Used to determine which answer keys to be pulling from.
 */

std::string getFileName(std::string file_path);

/**
 * @brief Checks to see if there is a directory at _dir_path_. If none exists, it creates on in place.
 * 
 * @pre The process running this function has permission to create a directory at _dir_path_ if need be.
 * 
 * @post The directory _dir_path_ exists at the specified location.
 */
void validate_directory(std::string dir_path);

/**
 * @brief Converts an array xin in sparse dense ordering to an array xout in
 * natural ordering
 */
void spdensetonatural(const double *xin,double *xout,int *idxn2sd_map,int nx);

/**
 * @brief Converts an array xin in natural ordering to an array xout in
 * sparse-dense ordering
 */
void naturaltospdense(const double *xin,double *xout,int *idxn2sd_map,int nx);

