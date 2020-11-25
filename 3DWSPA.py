from __future__ import division, print_function

import ntpath

import numpy.linalg as npl
import numpy as np
import sklearn
from scipy.sparse.csgraph import shortest_path
import sys
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import euclidean_distances


def contacts2distances(contacts):
    """ Infer distances from contact matrix"""
    dist_matrix = shortest_path(csgraph=contacts, method='J', directed=False)
    return dist_matrix


def distances2coordinates(distances):
    """ Infer coordinates from distances"""
    N = distances.shape[0]
    d_0 = []
    # pre-caching
    cache = {}
    for j in range(N):
        sumi = sum([distances[j, k] ** 2 for k in range(j + 1, N)])
        cache[j] = sumi

    # compute distances from center of mass
    sum2 = sum([cache[j] for j in range(N)])
    for i in range(N):
        sum1 = cache[i] + sum([distances[j, i] ** 2 for j in range(i + 1)])
        val = 1 / N * sum1 - 1 / N ** 2 * sum2
        d_0.append(val)

    # generate gram matrix
    gram = np.zeros(distances.shape)
    for row in range(distances.shape[0]):
        for col in range(distances.shape[1]):
            dists = d_0[row] ** 2 + d_0[col] ** 2 - distances[row, col] ** 2
            gram[row, col] = 1 / 2 * dists

    # extract coordinates from gram matrix
    coordinates = []
    vals, vecs = npl.eigh(gram)

    vals = vals[N - 3:]
    vecs = vecs.T[N - 3:]

    # must all be positive for PSD (positive semidefinite) matrix
    # same eigenvalues might be small -> exact embedding does not exist
    # fix by replacing all but largest 3 eigvals by 0
    # better if three largest eigvals are separated by large spectral gap

    for val, vec in zip(vals, vecs):
        coord = vec * np.sqrt(val)
        coordinates.append(coord)

    return np.array(coordinates).T


def calculateEuclideanDistace(distance_matrix):
    ecuclidean_distance_matrix = euclidean_distances(distance_matrix, distance_matrix)
    jdismatricF = distance_matrix.flatten()
    edFlatten_matrix = ecuclidean_distance_matrix.flatten()
    return jdismatricF, edFlatten_matrix


def apply_3DWSPA(contacts):
    """ Apply algorithm to data in given file"""
    distances = contacts2distances(contacts)
    coordinates = distances2coordinates(distances)
    return distances, coordinates


def writeLogFile(inputfile, allprintValues, outputFile):
    f = open(outputFile + "_" + "log.txt", "w")
    f.write("Input Filename: " + str(inputfile) + "\n")
    f.write("Convert factor: " + str(allprintValues[0]) + "\n")
    f.write("AVG Spearman correlation Dist vs. Reconstructed Dist: " + str(allprintValues[1]) + "\n")
    f.write("AVG Pearson correlation Dist vs. Reconstructed Dist: " + str(allprintValues[2]) + "\n")
    f.write("AVG RMSE: " + str(allprintValues[3]) + "\n")
    f.close()
    print("Input Filename: " + str(inputfile))
    print("Convert factor: " + str(allprintValues[0]))
    print("AVG Spearman correlation Dist vs. Reconstructed Dist: " + str(allprintValues[1]))
    print("AVG Pearson correlation Dist vs. Reconstructed Dist: " + str(allprintValues[2]))
    print("AVG RMSE: " + str(allprintValues[3]) + "\n")


def main():
    """ Main function"""
    # Read the input file and parameters
    global xyz
    my_dict = {}
    inputfile = sys.argv[1]
    f = open(inputfile, "r")
    for input_paramters in f:
        if input_paramters.startswith("#") or input_paramters == "\n":
            continue
        else:
            x, y = input_paramters.split('=')
            my_dict[x.strip()] = y.strip()

    filename = my_dict['INPUT_FILE']
    filen, form = ntpath.basename(filename).split('.')
    outputFilePath = my_dict['OUTPUT_FOLDER'] + filen + "_Output"
    conversion_factor_range = np.linspace(0, 2, 21)
    contacts = np.genfromtxt(filename, delimiter="\t")
    allSpearManCodd = []
    allprintValues = []
    # Convert the input N * N Matrix If's using conversion factor to distances
    for conversionFactor in conversion_factor_range:
        xyz = [[0 if item == 0 else 1 / pow(item, float(conversionFactor)) for item in subl] for subl in contacts]
        xyz = np.array(xyz)
        distance_matrixj, rec_coords = apply_3DWSPA(xyz)

        # Calculate the accuracy of the model using dSCC,dPSC,dRMSE
        # The value for dPCC and dSCC range from − 1 to + 1, where the higher values are preferred
        # Lower values of dRMSE is preferred
        flat1, flat2 = calculateEuclideanDistace(distance_matrixj)

        # distance Spearman Correlation Coefficient
        spearMancoff, pvalue = spearmanr(flat1, flat2)

        # distance Pearson correlation coefficient
        pearsonCoff, ppvalue = pearsonr(flat1, flat2)

        # distance root mean square error
        meanSqError = sklearn.metrics.mean_squared_error(flat1, flat2)
        drmse = np.math.sqrt(meanSqError)
        allSpearManCodd.append(spearMancoff)
        allprintValues.append([conversionFactor, spearMancoff, pearsonCoff, drmse, rec_coords])

    maxIndexValue = (allSpearManCodd.index(max(allSpearManCodd)))
    valuesForPdbFiles = allprintValues[maxIndexValue]
    # Write Log File
    writeLogFile(filename, allprintValues[maxIndexValue], outputFilePath)
    # Create the output pdb file
    make_pdb_file_00(outputFilePath, valuesForPdbFiles[4] * 200)


def make_pdb_file_00(outputfilepath, X):
    '''
        taks a 3xn array and converts it into a pdb format which can be read by pyMol.
        Referred from SIMBA3D
        '''
    with open(outputfilepath +".pdb", 'w') as result:
        for ii in range(len(X)):
            result.write('ATOM  ')
            result.write('{: 5d}'.format(ii + 1))
            result.write('   CA MET A' + str(ii + 1).ljust(8))
            result.write('{: 8.3f}'.format(X[ii, 0]))
            result.write('{: 8.3f}'.format(X[ii, 1]))
            result.write('{: 8.3f}'.format(X[ii, 2]))
            result.write('  0.20 10.00\n')
        for ii in range(len(X) - 1):
            result.write('CONECT')
            result.write('{: 5d}'.format(ii + 1))
            result.write('{: 5d}'.format(ii + 2) + '\n')
        result.write('END  ')


if __name__ == '__main__':
    main()
