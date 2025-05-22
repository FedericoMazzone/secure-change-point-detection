#include "diff-privacy.h"
#include "serial.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-ptxt.h"

#include "math/chebyshev.h"

#include <arpa/inet.h>
#include <fcntl.h>      // for O_RDWR
#include <filesystem>
#include <libgen.h>     // for dirname
#include <regex>
#include <sys/mman.h>   // for shm_open, mmap, munmap, MAP_SHARED, PROT_READ, PROT_WRITE, MAP_FAILED
#include <sys/stat.h>   // for S_IRUSR, S_IWUSR
#include <sys/wait.h>   // for waitpid()
#include <unistd.h>     // for fork(), execl()


using namespace lbcrypto;


// Generate coefficients of compare's polynomial approximation of given degree.
std::vector<double> generateCompareCoeffs
(
    double a,
    double b,
    uint32_t degree,
    double error = 0.00001
)
{
    return EvalChebyshevCoefficients(
        [error](double x) -> double { 
            if      (x > error)   return 1;
            else if (x >= -error) return 0.5;
            else                  return 0;
        },
        a, b, degree
    );
}


// Generate coefficients of indicator's polynomial approximation of degree k-1.
std::vector<double> generateIndicatorCoeffs
(
    int k,
    uint32_t degree
)
{
    if (degree == -1U)
    {
        // Resultant coefficients vector (initially size k, all zeros)
        std::vector<double> coeffs(k, 0.0);

        coeffs[0] = 1.0;

        // Iterate over each factor (x - i)
        for (int i = 1; i <= k - 1; i++)
        {
            for (int j = i; j >= 1; j--)
                coeffs[j] = coeffs[j - 1] - coeffs[j] * (i + 0.5);
            coeffs[0] = coeffs[0] * (-i - 0.5);
        }

        // Compute evaluation in 0.5
        double norm = coeffs[0];
        for (int i = 1; i <= k - 1; i++)
            norm += coeffs[i] * std::pow(0.5, i);
        
        // Normalize coefficients
        for (double& coeff : coeffs)
            coeff /= norm;

        return coeffs;
    }
    else
    {
        return EvalChebyshevCoefficients(
            [](double x) -> double { 
                if   (x > 1.0) return 0.0;
                else           return 1.0;
            },
            -0.5, k + 0.5, degree
        );
    }
}


std::vector<std::vector<std::vector<double>>> lloyd
(
    const std::vector<std::vector<double>> xA,
    const std::vector<std::vector<double>> xB,
    const size_t k,
    std::vector<std::vector<double>> cA,
    std::vector<std::vector<double>> cB,
    size_t &numRounds
)
{
    const size_t n = xA.size();
    const size_t ndimA = xA[0].size();
    const size_t ndimB = xB[0].size();
    const size_t ndim = ndimA + ndimB;
    bool hasConverged = false;
    std::vector<std::vector<double>> gA(k);
    std::vector<std::vector<double>> gB(k);
    std::vector<size_t> t(n);
    while (!hasConverged)
    {
        std::cout << "Centers Alice: " << cA << std::endl;
        std::cout << "Centers Bob  : " << cB << std::endl;

        // Centers to clusters
        for (size_t i = 0; i < n; i++)
        {
            size_t minIndex = 0;
            double minDist = ndim;
            for (size_t j = 0; j < k; j++)
            {
                double currDist = 0.0;
                for (size_t d = 0; d < ndimA; d++)
                    currDist += std::pow(cA[j][d] - xA[i][d], 2);
                for (size_t d = 0; d < ndimB; d++)
                    currDist += std::pow(cB[j][d] - xB[i][d], 2);
                if (currDist < minDist)
                {
                    minIndex = j;
                    minDist = currDist;
                }
            }
            t[i] = minIndex;
        }

        // Clusters to centers
        std::vector<size_t> clusterSize(k);
        for (size_t j = 0; j < k; j++)
        {
            gA[j] = std::vector<double>(ndimA, 0.0);
            gB[j] = std::vector<double>(ndimB, 0.0);
        }
        for (size_t i = 0; i < n; i++)
        {
            size_t clusterIndex = t[i];
            clusterSize[clusterIndex] += 1;
            for (size_t d = 0; d < ndimA; d++)
                gA[clusterIndex][d] += xA[i][d];
            for (size_t d = 0; d < ndimB; d++)
                gB[clusterIndex][d] += xB[i][d];
        }
        for (size_t j = 0; j < k; j++)
        {
            if (clusterSize[j] > 0)
            {
                for (size_t d = 0; d < ndimA; d++)
                    gA[j][d] /= clusterSize[j];
                for (size_t d = 0; d < ndimB; d++)
                    gB[j][d] /= clusterSize[j];
            }
        }

        hasConverged = true;
        for (size_t j = 0; j < k; j++)
        {
            if (gA[j] != cA[j] or gB[j] != cB[j])
                hasConverged = false;
            cA[j] = gA[j];
            cB[j] = gB[j];
        }
        
        numRounds++;
    }

    return {cA, cB};
}


// get the number of NUMA nodes
int getNumberOfNumaNodes()
{
    std::array<char, 128> buffer;
    std::string result;
    
    // Execute the "numactl --hardware" command
    std::unique_ptr<FILE, int(*)(FILE*)> pipe(popen("numactl --hardware", "r"), static_cast<int(*)(FILE*)>(pclose));
    
    if (!pipe)
    {
        std::cerr << "Failed to run command." << std::endl;
        return -1;
    }

    // Read the command output
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();

    // Use a regex to search for "available: X nodes"
    std::regex re("available:\\s*(\\d+)\\s*nodes");
    std::smatch match;
    
    if (std::regex_search(result, match, re) && match.size() > 1)
        return std::stoi(match.str(1));
    else
    {
        std::cerr << "Failed to find NUMA node information." << std::endl;
        return -1;
    }
}


std::string getExecutablePath()
{
    char result[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    if (count != -1)
    {
        result[count] = '\0';
        return std::string(result);
    }
    return "";
}


std::string getExecutableDir()
{
    std::string execPath = getExecutablePath();
    if (!execPath.empty())
    {
        char* execDir = dirname(&execPath[0]);
        return std::string(execDir);
    }
    return "";
}


// Wait for all child processes to be ready
void wait
(
    const size_t numCiphertexts,
    const std::string &storagePath,
    const std::string &label
)
{
    std::set<size_t> settingUpProcs; // set of child processes stil setting up
    for (size_t v = 0; v < numCiphertexts; v++)
        settingUpProcs.insert(v);
    while (settingUpProcs.size() > 0)
    {
        // Track processes that are done setting up
        for (auto it = settingUpProcs.begin(); it != settingUpProcs.end(); )
        {
            size_t v = *it;
            if (doesFileExist(storagePath + "/" + label + "-" + std::to_string(v)))
            {
                // Remove process from setting up processes set
                it = settingUpProcs.erase(it);

                deleteFile(storagePath + "/" + label + "-" + std::to_string(v));
            }
            else it++;
        }
    }
}


// Read and aggregate child process' output
void processOutput
(
    const std::string& storagePath,
    Ciphertext<DCRTPoly>& result,
    const size_t v,
    char* sharedMem,
    const size_t serializedSize
)
{
    // Retrieve serialized output from shared memory
    std::string result_s(sharedMem + serializedSize * v, serializedSize);

    // Deserialize the output and aggregate it with the others
    result += Serial::deserialize<Ciphertext<DCRTPoly>>(result_s);

    // Set child process' output as not ready for next iteration
    deleteFile(storagePath + "/output-status-" + std::to_string(v));
}


// Get the maximum CPU temperature using the 'sensors' command
float getMaxTemperature()
{
    // Run the 'sensors' command and capture the output
    FILE* pipe = popen("sensors | grep 'Core' | awk '{print $3}' | sed 's/+//;s/°C//'", "r");
    if (!pipe)
    {
        std::cerr << "Failed to run sensors command." << std::endl;
        return -1;
    }
    char buffer[128];
    std::string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr)
        result += buffer;
    pclose(pipe);

    // Parse the temperatures and find the maximum
    std::istringstream ss(result);
    float temp, maxTemp = -1000;
    while (ss >> temp)
        if (temp > maxTemp)
            maxTemp = temp;
    return maxTemp;
}


// Wait until the CPU temperature is below a specified threshold
void waitForOptimalTemperature
(
    float threshold
)
{
    while (true)
    {
        float maxTemp = getMaxTemperature();
        std::cout << "Current max temperature: " << maxTemp << "°C" << std::endl;
        if (maxTemp < threshold)
        {
            std::cout << "Temperature is below threshold (" << threshold << "°C). Proceeding..." << std::endl;
            break;
        }
        else
        {
            std::cout << "Waiting for the temperature to cool down..." << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(2));
        }
    }
}


std::vector<uint> parseCommaSeparatedValues
(
    const std::string& input
)
{
    std::vector<uint> values;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ','))
        try
        {
            values.push_back(std::stoi(token));
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "Invalid integer in input: " << token << std::endl;
        }

    return values;
}


// Normalizing points in [0,1]^numDim
void normalize
(
    std::vector<std::vector<double>> &x
)
{
    for (size_t d = 0; d < x[0].size(); d++)
    {
        double xmin = x[0][d];
        double xmax = x[0][d];
        for (size_t i = 1; i < x.size(); i++)
        {
            if (x[i][d] < xmin)
                xmin = x[i][d];
            if (x[i][d] > xmax)
                xmax = x[i][d];
        }
        for (size_t i = 0; i < x.size(); i++)
            x[i][d] = (x[i][d] - xmin) / (xmax - xmin);
    }
}


// Hungarian algorithm
// Andrew Lopatin's implementation
// http://e-maxx.ru/algo/assignment_hungary#6
int hungarian
(
    std::vector<std::vector<int>>& costMatrix,
    std::vector<int>& assignment
)
{
    int n = costMatrix.size();
    for (int i = 0; i < n; i++)
        costMatrix[i].insert(costMatrix[i].begin(), 0);
    costMatrix.insert(costMatrix.begin(), std::vector<int>(n + 1, 0));

    std::vector<int> u(n + 1), v(n + 1), p(n + 1), way(n + 1);
    for (int i = 1; i <= n; i++)
    {
        p[0] = i;
        int j0 = 0;
        std::vector<int> minv(n + 1, std::numeric_limits<int>::max());
        std::vector<char> used(n + 1, false);
        do
        {
            used[j0] = true;
            int i0 = p[j0], delta = std::numeric_limits<int>::max(), j1 = -1;
            for (int j = 1; j <= n; j++)
                if (!used[j])
                {
                    int cur = costMatrix[i0][j] - u[i0] - v[j];
                    if (cur < minv[j])
                        minv[j] = cur, way[j] = j0;
                    if (minv[j] < delta)
                        delta = minv[j], j1 = j;
                }
            for (int j = 0; j <= n; j++)
                if (used[j])
                    u[p[j]] += delta, v[j] -= delta;
                else
                    minv[j] -= delta;
            j0 = j1;
        }
        while (p[j0] != 0);
        do
        {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        }
        while (j0);
    }

    // Generate result from dual variables
    for (int j = 1; j <= n; j++)
        assignment[p[j]] = j;

    // Calculate the optimal cost
    int optimalCost = -v[0];
    return optimalCost;
}


double computeLoss
(
    const std::vector<std::vector<double>> &xA,
    const std::vector<std::vector<double>> &xB,
    const std::vector<std::vector<double>> &cA,
    const std::vector<std::vector<double>> &cB
)
{
    const size_t n = xA.size();
    const size_t k = cA.size();
    const size_t ndimA = xA[0].size();
    const size_t ndimB = xB[0].size();
    const size_t ndim = ndimA + ndimB;
    double loss = 0.0;
    
    for (size_t i = 0; i < n; i++)
    {
        double minDist = ndim;
        for (size_t j = 0; j < k; j++)
        {
            double currDist = 0.0;
            for (size_t d = 0; d < ndimA; d++)
                currDist += std::pow(cA[j][d] - xA[i][d], 2);
            for (size_t d = 0; d < ndimB; d++)
                currDist += std::pow(cB[j][d] - xB[i][d], 2);
            if (currDist < minDist)
                minDist = currDist;
        }
        loss += minDist;
    }

    return loss / n;
}


// compute clusters from centroids
std::vector<int> computeClusters
(
    const std::vector<std::vector<double>> &xA,
    const std::vector<std::vector<double>> &xB,
    const std::vector<std::vector<double>> &cA,
    const std::vector<std::vector<double>> &cB
)
{
    const size_t n = xA.size();
    const size_t k = cA.size();
    const size_t ndimA = xA[0].size();
    const size_t ndimB = xB[0].size();
    const size_t ndim = ndimA + ndimB;
    std::vector<int> clusters(n);

    for (size_t i = 0; i < n; i++)
    {
        size_t minIndex = 0;
        double minDist = ndim;
        for (size_t j = 0; j < k; j++)
        {
            double currDist = 0.0;
            for (size_t d = 0; d < ndimA; d++)
                currDist += std::pow(cA[j][d] - xA[i][d], 2);
            for (size_t d = 0; d < ndimB; d++)
                currDist += std::pow(cB[j][d] - xB[i][d], 2);
            if (currDist < minDist)
            {
                minIndex = j;
                minDist = currDist;
            }
        }
        clusters[i] = minIndex;
    }

    return clusters;
}


double clusterAccuracy
(
    const std::vector<std::vector<double>> &xA,
    const std::vector<std::vector<double>> &xB,
    const std::vector<std::vector<double>> &cA,
    const std::vector<std::vector<double>> &cB,
    const std::vector<int> &groundTruth
)
{
    const size_t n = xA.size();
    const size_t k = cA.size();
    const size_t ndimA = xA[0].size();
    const size_t ndimB = xB[0].size();
    const size_t ndim = ndimA + ndimB;

    // Building costMatrix, where costMatrix[i][j] is the number of
    // points with ground-truth cluster i and predicted cluster j
    std::vector<std::vector<int>> costMatrix(k, std::vector<int>(k, 0));
    for (size_t i = 0; i < n; i++)
    {
        size_t minIndex = 0;
        double minDist = ndim;
        for (size_t j = 0; j < k; j++)
        {
            double currDist = 0.0;
            for (size_t d = 0; d < ndimA; d++)
                currDist += std::pow(cA[j][d] - xA[i][d], 2);
            for (size_t d = 0; d < ndimB; d++)
                currDist += std::pow(cB[j][d] - xB[i][d], 2);
            if (currDist < minDist)
            {
                minIndex = j;
                minDist = currDist;
            }
        }
        costMatrix[groundTruth[i]][minIndex] -= 1;
    }

    // Apply Hungarian algorithm to find optimal assignment
    std::vector<int> assignment(k + 1);
    int optimalCost = hungarian(costMatrix, assignment);
    size_t correctlyClassified = -optimalCost;
    double accuracy = 1.0 * correctlyClassified / n;

    // Rearrange cost matrix
    std::vector<std::vector<int>> confusionMatrix(k, std::vector<int>(k, 0));
    for (size_t i = 1; i < assignment.size(); i++) {
        int col = assignment[i];  // Optimal column for row i
        for (size_t row = 1; row < costMatrix.size(); row++)
            confusionMatrix[row - 1][i - 1] = -costMatrix[row][col];  // Move the column
    }
    std::cout << "Cluster accuracy: " << accuracy << std::endl;
    std::cout << "Confusion matrix: " << confusionMatrix << std::endl;

    return accuracy;
}


// Rank a vector, handling duplicates by assigning the same rank
std::vector<int> rankVector
(
    const std::vector<int>& input
)
{
    // Create a vector of pairs where each pair contains the value and its index
    std::vector<std::pair<int, int>> valueIndexPairs;
    for (size_t i = 0; i < input.size(); i++)
        valueIndexPairs.emplace_back(input[i], i);

    // Sort the vector of pairs by the value (first element of the pair)
    std::sort(valueIndexPairs.begin(), valueIndexPairs.end());

    // Create a result vector initialized to the same size as the input
    std::vector<int> result(input.size());

    // Assign ranks to the elements, handling duplicates
    int rank = 0;
    result[valueIndexPairs[0].second] = rank;
    for (size_t i = 1; i < valueIndexPairs.size(); i++)
    {
        // If the current value is the same as the previous, assign the same rank
        if (valueIndexPairs[i].first != valueIndexPairs[i - 1].first)
            rank++;
        result[valueIndexPairs[i].second] = rank;
    }

    return result;
}


/**
 * @brief Load a dataset from a file.
 *
 * Read a dataset from the specified CSV file, separate each data point into
 * two subsets of features based on specified dimension indices, and optionally
 * extract the ground-truth labels (as last column).
 *
 * @param filepath Path to the dataset file. Each line in the file represents a
 *        data point, with features separated by commas.
 * @param dimsA Vector specifying indices of features to extract into subset A.
 * @param dimsB Vector specifying indices of features to extract into subset B.
 * @param testClusters Boolean flag indicating if the last value in each line
 *        represents a ground-truth label. If true, the function will store
 *        these labels in `groundTruth`.
 * @param[out] n The number of data points loaded from the file.
 * @param[out] xA 2D vector to store features specified by `dimsA` for each
 *             data point.
 * @param[out] xB 2D vector to store features specified by `dimsB` for each
 *             data point.
 * @param[out] groundTruth Vector to store ground-truth cluster labels if
 *             `testClusters` is true. The labels are normalized to the range
 *             [0, k-1].
 *
 * @throws std::runtime_error If the file cannot be opened.
 *
 * @note The data points in `xA` and `xB` are normalized to fit within the
 *       [0,1] range. This function assumes that the file format is CSV with
 *       numeric values, and that each row has enough values to accommodate all
 *       indices specified in `dimsA` and `dimsB`.
 */
void loadDataset(
    const std::string &filepath,
    const std::vector<uint> &dimsA,
    const std::vector<uint> &dimsB,
    const bool testClusters,
    size_t &n,
    std::vector<std::vector<double>> &xA,
    std::vector<std::vector<double>> &xB,
    std::vector<int> &groundTruth
)
{
    // Opening datapoints file
    std::ifstream infile(filepath);
    if (!infile)
        throw std::runtime_error("Error: Unable to open file " + filepath);

    // Reading points from file
    std::string line;
    while (std::getline(infile, line))
    {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;

        // Split each line by comma
        while (std::getline(ss, value, ','))
            row.push_back(std::stod(value));

        // Append the values to the corresponding vectors
        std::vector<double> rowA;
        for (size_t dim : dimsA)
            rowA.push_back(row[dim]);
        xA.push_back(rowA);
        std::vector<double> rowB;
        for (size_t dim : dimsB)
            rowB.push_back(row[dim]);
        xB.push_back(rowB);

        // Store ground-truth label if available
        if (testClusters)
            groundTruth.push_back((int) row[row.size() - 1]);

        n++;
    }
    infile.close();

    // Shift values in groundTruth to be in 0..k-1
    if (testClusters)
        groundTruth = rankVector(groundTruth);

    // Normalize data in [0,1]^ndim
    normalize(xA);
    normalize(xB);
}


void initializeCentroidsRandom
(
    const size_t k,
    const size_t ndimA,
    const size_t ndimB,
    std::vector<std::vector<double>> &cA,
    std::vector<std::vector<double>> &cB
)
{
    srand((unsigned)time(NULL));
    for (size_t j = 0; j < k; j++)
    {
        for (size_t d = 0; d < ndimA; d++)
            cA[j].push_back((double) rand() / (double) RAND_MAX);
        for (size_t d = 0; d < ndimB; d++)
            cB[j].push_back((double) rand() / (double) RAND_MAX);
    }
}


// for testing purposes only
void initializeCentroidsGoodStart
(
    const size_t k,
    const std::vector<int> &groundTruth,
    const std::vector<std::vector<double>> &xA,
    const std::vector<std::vector<double>> &xB,
    std::vector<std::vector<double>> &cA,
    std::vector<std::vector<double>> &cB
)
{
    for (size_t j = 0; j < k; j++)
    {
        auto it = std::find(groundTruth.begin(), groundTruth.end(), j);
        size_t i = std::distance(groundTruth.begin(), it);
        cA[j] = xA[i];
        cB[j] = xB[i];
    }
}


void initializeCentroidsGrid
(
    const size_t k,
    const size_t ndimA,
    const size_t ndimB,
    std::vector<std::vector<double>> &cA,
    std::vector<std::vector<double>> &cB
)
{
    const size_t ndim = ndimA + ndimB;

    // Determine grid size per dimension (approximate nth root of k)
    size_t gridSize = static_cast<size_t>(std::ceil(std::pow(k, 1.0 / ndim)));

    // Calculate the total number of grid points
    size_t totalGridPoints = static_cast<size_t>(std::pow(gridSize, ndim));

    // Randomly select k unique indices from the total grid points
    std::unordered_set<size_t> uniqueIndices;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, totalGridPoints - 1);
    while (uniqueIndices.size() < k)
        uniqueIndices.insert(dist(gen));
    
    // Convert each selected index into n-dimensional coordinates
    double step = 1.0 / gridSize;
    size_t j = 0;
    for (size_t index : uniqueIndices)
    {
        // Convert a linear index to n-dimensional grid coordinates
        std::vector<size_t> gridCoords(ndim);
        for (size_t d = 0; d < ndim; d++)
        {
            gridCoords[d] = index % gridSize;
            index /= gridSize;
        }

        // Map grid coordinates to [0, 1] range
        cA[j] = std::vector<double>(ndimA);
        for (size_t d = 0; d < ndimA; d++)
            cA[j][d] = (gridCoords[d] + 0.5) * step;
        cB[j] = std::vector<double>(ndimB);
        for (size_t d = 0; d < ndimB; d++)
            cB[j][d] = (gridCoords[d + ndimA] + 0.5) * step;
        j++;
    }
}


int createServer
(
    const char *addr,
    const uint16_t port,
    int &server_fd
)
{
    // Create socket
    server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (server_fd == -1)
        throw std::runtime_error("Server Creation: Socket creation failed!");

    // Setup address structure
    sockaddr_in address;
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = inet_addr(addr);
    address.sin_port = htons(port);

    // Bind the socket
    if (bind(server_fd, (struct sockaddr *)&address, sizeof(address)) < 0)
    {
        close(server_fd);
        throw std::runtime_error("Server Creation: Bind failed!");
    }

    // Listen for incoming connections
    if (listen(server_fd, 3) < 0)
    {
        close(server_fd);
        throw std::runtime_error("Server Creation: Listen failed!");
    }

    std::cout << "Server is listening on " << addr << ":" << port << std::endl;

    // Accept a connection
    int clientSocket = accept(server_fd, nullptr, nullptr);
    if (clientSocket < 0)
    {
        close(server_fd);
        throw std::runtime_error("Server Creation: Accept failed!");
    }

    // Retrieve client information
    sockaddr_in clientAddress;
    socklen_t clientAddressLen = sizeof(clientAddress);
    if (getpeername(clientSocket, (struct sockaddr*)&clientAddress, &clientAddressLen) == 0)
    {
        // Convert client address and port to readable format
        char clientAddr[INET_ADDRSTRLEN];
        inet_ntop(AF_INET, &clientAddress.sin_addr, clientAddr, sizeof(clientAddr));
        int clientPort = ntohs(clientAddress.sin_port);

        std::cout << "Connected to client at " << clientAddr << ":" << clientPort << std::endl;
    }
    else
        std::cerr << "Failed to get client address!" << std::endl;
    
    return clientSocket;
}


int connectToServer
(
    const char *serverAddr,
    const uint16_t serverPort
)
{
    // Create socket
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == -1)
        throw std::runtime_error("Connection to Server: Socket creation failed!");

    // Setup server address
    sockaddr_in serverAddress;
    serverAddress.sin_family = AF_INET;
    serverAddress.sin_port = htons(serverPort);
    serverAddress.sin_addr.s_addr = inet_addr(serverAddr); // Loopback interface

    // Connect to the server
    if (connect(sock, (struct sockaddr *)&serverAddress, sizeof(serverAddress)) < 0)
    {
        close(sock);
        throw std::runtime_error("Connection to Server: Connection failed!");
    }

    return sock;
}


int commSize = 0;

void sendMessage
(
    const int socket,
    const std::string &message
)
{
    // Get message length
    uint32_t length = message.size();
    uint32_t lengthNetwork = htonl(length); // Convert to network byte order
    commSize += length + 4;

    // Send length header
    if (send(socket, &lengthNetwork, sizeof(lengthNetwork), 0) < 0)
        throw std::runtime_error("Issue with sending message length!");

    // Send message
    if (send(socket, message.c_str(), length, 0) < 0)
        throw std::runtime_error("Issue with sending message!");
}


std::string receiveMessage
(
    const int socket
)
{
    // Read length header
    uint32_t lengthNetwork;
    if (recv(socket, &lengthNetwork, sizeof(lengthNetwork), 0) <= 0)
        throw std::runtime_error("Issue with receiving message length!");

    // Convert length from network byte order to host byte order
    uint32_t length = ntohl(lengthNetwork);
    commSize += length + 4;

    // Read message data
    std::vector<char> buffer(length);
    uint32_t bytesRead = 0;
    while (bytesRead < length)
    {
        int result = recv(socket, &buffer[bytesRead], length - bytesRead, 0);
        if (result <= 0)
            throw std::runtime_error("Issue with receiving message!");
        bytesRead += result;
    }

    return std::string(buffer.begin(), buffer.end());
}


void writeToBinFile
(
    const std::string &filename,
    const std::string &content
)
{
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile)
        throw std::runtime_error("Error opening file for writing: " + filename);
    outFile << content;
    outFile.close();
}


int main(int argc, char* argv[])
{

    bool verbose = false;

    // File with datapoints
    std::string filepath = "";

    // Number of clusters
    size_t k = 0;

    // Number of rounds
    size_t numRounds = -1UL;

    // Differential Privacy parameters
    double sensitivity    = 1.0;  // sensitivity of the new centers
    double epsilonTotal   = 1.0;  // epsilon (privacy parameter)
    double deltaTotal     = 1e-6; // delta (failure probability)

    // Convergence parameters
    const double convergenceParam = 1.0;
    const size_t convergenceWindowSize = 2;

    // Approximation degrees
    usint compareDepth = 10;
    usint indicatorDepth = -1U;
    // k = 2 -> no indicator used, use default -1U
    // k = 3 -> 4 (standard is better, use default -1U)
    // k = 5 -> 4 (standard is better, use default -1U)
    // k = 8 -> 4 (standard is better, use default -1U)
    // k = 15 -> 5 (standard does not work, use chebyshev)

    // Define the maximum acceptable CPU temperature to proceed (in °C)
    const float temperatureThreshold = 50.0;

    // Dimension distribution
    std::vector<uint> dimsA = {0};
    std::vector<uint> dimsB = {1};

    // Are ground-truth cluster ids available
    bool testClusters = false;

    // Party ID for simulating communication on loopback interface
    int partyId = -1;


    // Parse args
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        
        if (arg == "-v" || arg == "--verbose")
            verbose = true;
        else if (arg == "-f" || arg == "--file")
            if (i + 1 < argc)
                filepath = argv[++i];
            else
            {
                std::cerr << "Error: --file requires a filepath argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-k" || arg == "--clusters")
            if (i + 1 < argc)
                k = std::stoul(argv[++i]);
            else
            {
                std::cerr << "Error: --clusters requires a num_clusters argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-cd" || arg == "--cmp-degree")
            if (i + 1 < argc)
                compareDepth = std::stoi(argv[++i]);
            else
            {
                std::cerr << "Error: --cmp-degree requires a degree argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-id" || arg == "--ind-degree")
            if (i + 1 < argc)
                indicatorDepth = std::stoi(argv[++i]);
            else
            {
                std::cerr << "Error: --ind-degree requires a degree argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-r" || arg == "--rounds")
            if (i + 1 < argc)
                numRounds = std::stoul(argv[++i]);
            else
            {
                std::cerr << "Error: --rounds requires a num_rounds argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-e" || arg == "--epsilon")
            if (i + 1 < argc)
                epsilonTotal = std::stod(argv[++i]);
            else
            {
                std::cerr << "Error: --epsilon requires a privacy budget argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-dimsA")
            if (i + 1 < argc)
                dimsA = parseCommaSeparatedValues(argv[++i]);
            else
            {
                std::cerr << "Error: -dimsA requires a dimensions argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-dimsB")
            if (i + 1 < argc)
                dimsB = parseCommaSeparatedValues(argv[++i]);
            else
            {
                std::cerr << "Error: -dimsB requires a dimensions argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-t" || arg == "--test-clusters")
            testClusters = true;
        else if (arg == "-p" || arg == "--party")
            if (i + 1 < argc)
            {
                partyId = std::stoi(argv[++i]);
                if (partyId != 0 && partyId != 1)
                {
                    std::cerr << "Error: party id can only be 0 or 1."
                              << std::endl;
                    return 1;
                }
            }
            else
            {
                std::cerr << "Error: --party requires a party id argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-h" || arg == "--help")
        {
            std::cout << "Usage: clustering "
                      << "-f|--file <filepath> "
                      << "-k|--clusters <num_clusters> "
                      << "[-r|--rounds <num_rounds>] "
                      << "[-dimsA <dim1,dim2,...>] "
                      << "[-dimsB <dim1,dim2,...>] "
                      << "[-cd|--cmp-degree <degree>] "
                      << "[-id|--ind-degree <degree>] "
                      << "[-e|--epsilon <privacy_budget>] "
                      << "[-t|--test-clusters] "
                      << "[-v|--verbose] "
                      << "[-h|--help] "
                      << "[-p|--party <party_id>]"
                      << std::endl;
            return 0;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }

    // Checking for required args
    if (k == 0)
    {
        std::cerr << "Missing number of cluster k." << std::endl;
        return 1;
    }
    if (filepath.length() == 0)
    {
        std::cerr << "Missing datapoints filepath." << std::endl;
        return 1;
    }

    size_t n = 0;                           // Number of points
    std::vector<std::vector<double>> xA;    // Features of Alice
    std::vector<std::vector<double>> xB;    // Features of Bob
    std::vector<int> groundTruth;           // Cluster groun-truth labels

    loadDataset(
        filepath, dimsA, dimsB, testClusters,
        n, xA, xB, groundTruth
    );

    deltaTotal = 1.0 / n;

    const int numNumaNodes = getNumberOfNumaNodes();
    const std::string subprocPath = getExecutableDir() + "/clustering-sub" + ((k == 2) ? "-k2" : "");

    const bool indCheby = indicatorDepth != -1U;

    std::cout << std::endl <<
        "K-Means clustering via Lloyd algorithm on vertically partitioned data"
        << std::endl << std::endl;
    std::cout << "n : " << n << std::endl;
    std::cout << "k : " << k << std::endl << std::endl;
    std::cout << "Dimensions Alice: " << dimsA << std::endl;
    std::cout << "Dimensions Bob  : " << dimsB << std::endl << std::endl;
    // std::cout << "Datapoints Alice: " << xA << std::endl;
    // std::cout << "Datapoints Bob  : " << xB << std::endl << std::endl;
    std::cout << "Comparison depth: " << compareDepth << std::endl;
    std::cout << "Indicator depth : " << (indCheby ? std::to_string(indicatorDepth) : "equispaced") << std::endl << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;
    std::cout << "Number of NUMA nodes: " << numNumaNodes << std::endl;
    std::cout << "Expected sub-process executable path: " << subprocPath << std::endl;
    std::cout << std::endl;

    if (partyId != -1)
        std::cout << "Running as party " << partyId << std::endl << std::endl;
    
    // Initialize timers
    auto t = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSeconds;
    std::cout << std::fixed << std::setprecision(3);
    #undef TOC
    #define TOC(t) (elapsedSeconds = timeNow() - t).count()

    // Initialize server/client
    const char* address = "127.0.0.1";
    const uint16_t port = 1212;
    int server_fd = 0;
    int socket = 0;
    if (partyId == 0)
        socket = createServer(address, port, server_fd);
    else if (partyId == 1)
        socket = connectToServer(address, port);
    
    if (partyId != 1)
        TIC(t);

    
    ////// INITIALIZE CENTROIDS

    const size_t ndimA = dimsA.size();
    const size_t ndimB = dimsB.size();
    const size_t ndim = ndimA + ndimB;
    std::vector<std::vector<double>> cA(k), cB(k);

    // initializeCentroidsRandom(k, ndimA, ndimB, cA, cB);
    // initializeCentroidsGoodStart(k, groundTruth, xA, xB, cA, cB);
    initializeCentroidsGrid(k, ndimA, ndimB, cA, cB);



    ////////////////////////////////////////////////////////////////////////
    //                   LLoyd's Algorithm in Plaintext                   //
    ////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<std::vector<double>>> expectedCenters;
    if (partyId != 1)
    {
        std::cout << "Running LLoyd's Algorithm in Plaintext" << std::endl;
        size_t numRoundsPtxt = 0;
        expectedCenters = lloyd(xA, xB, k, cA, cB, numRoundsPtxt);
        std::cout << std::endl;
    }



    ////////////////////////////////////////////////////////////////////////
    //                     Differential Privacy Setup                     //
    ////////////////////////////////////////////////////////////////////////

    double epsilonRound = 0.0;
    double epsilonRoundAdv = 0.0;
    double deltaRound = 0.0;
    bool usingAdvComp = false;
    double noiseScale = 0.0;
    if (partyId != 1)
    {
        epsilonRound = DP::computePerRoundEpsilon(epsilonTotal, numRounds);
        epsilonRoundAdv = DP::computePerRoundEpsilonAdv(epsilonTotal, deltaTotal, numRounds);
        deltaRound = DP::computePerRoundDelta(deltaTotal, numRounds);
        usingAdvComp = epsilonRoundAdv > epsilonRound;
        noiseScale = DP::calculateGaussianNoiseScale(usingAdvComp ? epsilonRoundAdv : epsilonRound, deltaRound, sensitivity);
        std::cout << "DIFFERENTIAL PRIVACY PARAMETERS"                                << std::endl;
        std::cout << "Upper bound on the number of rounds      : " << numRounds       << std::endl;
        std::cout << "Sensitivity of the new centers           : " << sensitivity     << std::endl;
        std::cout << "Total privacy budget (epsilon)           : " << epsilonTotal    << std::endl;
        std::cout << "Total privacy loss (delta)               : " << std::scientific << deltaTotal << std::fixed << std::endl;
        std::cout << "Per-round epsilon (simple composition)   : " << epsilonRound    << std::endl;
        std::cout << "Per-round epsilon (advanced composition) : " << epsilonRoundAdv << std::endl;
        std::cout << "Using " << (usingAdvComp ? "advanced composition." : "simple composition.") << std::endl;
        std::cout << "Per-round delta (simple composition)     : " << std::scientific << deltaRound << std::fixed << std::endl;
        std::cout << "Gaussian noise scale                     : " << noiseScale      << std::endl;
        std::cout << std::endl;
    }



    ////////////////////////////////////////////////////////////////////////
    //                  Computing Convergence Condition                   //
    ////////////////////////////////////////////////////////////////////////

    double convergenceRadius = 0.0;
    if (partyId != 1)
    {
        // max between differential privacy and polynomial approximation error
        convergenceRadius = convergenceParam * MAX(noiseScale, 30.0);
        std::cout << "CONVERGENCE CONDITION PARAMETERS"                    << std::endl;
        std::cout << "Convergence Parameter   : " << convergenceParam      << std::endl;
        std::cout << "Convergence Radius      : " << convergenceRadius     << " / number of cluster points" << std::endl;
        std::cout << "Convergence Window Size : " << convergenceWindowSize << std::endl;
        std::cout << std::endl;
    }



    ////////////////////////////////////////////////////////////////////////
    //                       Setting up CKKS scheme                       //
    ////////////////////////////////////////////////////////////////////////

    CryptoContext<DCRTPoly> cryptoContext;
    if (partyId != 0)
    {
        usint indicatorCost = 0;
        if (indicatorDepth == -1U)
        {
            if (k == 2)
                indicatorCost = 0;
            else if (k == 3)
                indicatorCost = 2;
            else if (k <= 6)
                indicatorCost = 3;
            else if (k <= 10)
                indicatorCost = 4;
            else if (k <= 15)
                indicatorCost = 5;
            else
            {
                std::cerr << "Indicator cost not computed for k " << k << std::endl;
                return 1;
            }
        }
        else
            indicatorCost = indicatorDepth;
        
        const usint integralPrecision   = std::floor(std::log2(n)) + 1;
        const usint decimalPrecision    = 32; // 60 - (int) integralPrecision;
        const usint multiplicativeDepth = 3 + compareDepth + indicatorCost + ((ndim > 2 && k > 2) ? 1 : 0);

        cryptoContext = generateCryptoContext(
            integralPrecision,
            decimalPrecision,
            multiplicativeDepth,
            0, false, 0, true
        );
    }

    std::string cryptoContextSer;
    if (partyId == 1)
        sendMessage(socket, Serial::serialize(cryptoContext));
    else if (partyId == 0)
    {
        cryptoContextSer = receiveMessage(socket);
        cryptoContext = Serial::deserialize<CryptoContext<DCRTPoly>>(cryptoContextSer);
    }

    const usint numSlots = cryptoContext->GetCryptoParameters()->GetEncodingParams()->GetBatchSize();
    const usint maxNumEl = (k == 2) ? MIN(numSlots, n) : MIN(numSlots / (k * k), n);
    const usint numCiphertexts = std::ceil(1.0 * n / maxNumEl);
    const usint maxNumElCompressed = (k == 2) ? numSlots : (numSlots / (k * k)) * (k * k);
    const usint numCompressedCiphertexts = std::ceil(1.0 * n / maxNumElCompressed);
    std::cout << "Maximum number of elements per ciphertext            : " << maxNumEl                  << std::endl;
    std::cout << "Number of ciphertexts used                           : " << numCiphertexts            << std::endl;
    std::cout << "Maximum number of elements per compressed ciphertext : " << maxNumElCompressed        << std::endl;
    std::cout << "Number of compressed ciphertexts used                : " << numCompressedCiphertexts  << std::endl;
    std::cout << std::endl;

    // Generating public/private key pair, relinearization, and rotation keys
    // << positive
    // >> negative
    // assuming k <= n
    PublicKey<DCRTPoly> publicKey;
    PrivateKey<DCRTPoly> secretKey;
    if (partyId != 0)
    {
        std::vector<int32_t> indices;
        if (k == 2)
        {
            for (size_t l = 0; l < LOG2(maxNumEl); l++)
                indices.push_back(1 << l);
            for (size_t d = 1; d < 2 * ndim + 1; d++)
                indices.push_back(-d);
        }
        else
        {
            for (int32_t j = 1; j < (int32_t) k; j++)
            {
                indices.push_back(j);
                indices.push_back(j * maxNumEl * k);
                indices.push_back(-j);
                indices.push_back(-(j * maxNumEl * k));
            }
            for (int32_t l = 0; l < (int32_t) std::floor(std::log2(maxNumEl)); l++)
            {
                indices.push_back((1 << l) * k);
                if ((maxNumEl >> l) & 1)
                    indices.push_back((maxNumEl - (maxNumEl & ((1 << (l + 1)) - 1))) * k);
            }
            for (size_t d = 1; d < MAX(ndimA, ndimB); d++)
                indices.push_back(-d * k);
        }
        KeyPair<DCRTPoly> keyPair = keyGeneration(
            cryptoContext,
            indices,
            numSlots,
            false, true
        );
        publicKey = keyPair.publicKey;
        secretKey = keyPair.secretKey;
    }

    std::string publicKeySer;
    std::string relinKeySer;
    std::string rotationKeysSer;
    if (partyId == 1)
    {
        sendMessage(socket, Serial::serialize(publicKey));
        sendMessage(socket, Serial::serialize(cryptoContext->GetAllEvalMultKeys().begin()->second[0]));
        sendMessage(socket, Serial::serialize(cryptoContext->GetAllEvalAutomorphismKeys().begin()->second));
    }
    else if (partyId == 0)
    {
        publicKeySer = receiveMessage(socket);
        relinKeySer = receiveMessage(socket);
        rotationKeysSer = receiveMessage(socket);
        publicKey = Serial::deserialize<PublicKey<DCRTPoly>>(publicKeySer);
        // auto relinKey = Serial::deserialize<EvalKey<DCRTPoly>>(receiveMessage(socket));
        // auto rotationKeys = Serial::deserialize<std::shared_ptr<std::map<usint, EvalKey<DCRTPoly>>>>(receiveMessage(socket));
        // cryptoContext->InsertEvalMultKey({relinKey});
        // cryptoContext->InsertEvalAutomorphismKey(rotationKeys);
    }
    else if (partyId == -1)
    {
        cryptoContextSer = Serial::serialize(cryptoContext);
        publicKeySer = Serial::serialize(publicKey);
        relinKeySer = Serial::serialize(cryptoContext->GetAllEvalMultKeys().begin()->second[0]);
        rotationKeysSer = Serial::serialize(cryptoContext->GetAllEvalAutomorphismKeys().begin()->second);
    }



    ////////////////////////////////////////////////////////////////////////
    //                Pre-Computations (data-independent)                 //
    ////////////////////////////////////////////////////////////////////////

    std::vector<double> cmpCoeffs;
    std::vector<double> indCoeffs;
    Ciphertext<DCRTPoly> zeroC;
    Plaintext plaintext;
    std::vector<double> plaintextVec;
    if (partyId != 1)
    {
        cmpCoeffs = generateCompareCoeffs(-1.0, 1.0, depth2degree(compareDepth));
        indCoeffs = generateIndicatorCoeffs(k, depth2degree(indicatorDepth));
        zeroC = cryptoContext->Encrypt(
            publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(std::vector<double> {0.0})
        );
    }



    ////////////////////////////////////////////////////////////////////////
    //                      Prepare Child Processes                       //
    ////////////////////////////////////////////////////////////////////////


    std::string storagePath = "data/tmp";
    size_t serializedSize = 0;
    size_t sharedMemorySize = 0;
    char* sharedMem;
    std::map<pid_t, size_t> pidToId; // Map PID to the corresponding ciphertext id
    if (partyId != 1)
    {
        ////// STORING DATA FOR CHILD PROCESSES

        std::cout << "Storing data for subprocesses... " << std::flush;
        if (!std::filesystem::exists(storagePath))
            std::filesystem::create_directory(storagePath);

        Serial::serializeToFile(n, storagePath + "/n");
        Serial::serializeToFile(k, storagePath + "/k");
        writeToBinFile(storagePath + "/crypto-context", cryptoContextSer);
        // Serial::serializeToFile(cryptoContext, storagePath + "/crypto-context");
        writeToBinFile(storagePath + "/public-key", publicKeySer);
        // Serial::serializeToFile(publicKey, storagePath + "/public-key");
        // Serial::serializeToFile(secretKey, storagePath + "/secret-key");
        writeToBinFile(storagePath + "/relin-key", relinKeySer);
        writeToBinFile(storagePath + "/rotation-keys", rotationKeysSer);
        // Serial::serializeToFile(cryptoContext->GetAllEvalMultKeys().begin()->second[0], storagePath + "/relin-key");
        // Serial::serializeToFile(cryptoContext->GetAllEvalAutomorphismKeys().begin()->second, storagePath + "/rotation-keys");
        Serial::serializeToFile(cmpCoeffs, storagePath + "/cmp-coeffs");
        Serial::serializeToFile(indCoeffs, storagePath + "/ind-coeffs");
        deleteFile(storagePath + "/start-preparation");
        deleteFile(storagePath + "/start-computing");
        deleteFile(storagePath + "/stop");
        for (size_t v = 0; v < numCiphertexts; v++)
        {
            deleteFile(storagePath + "/proc-ready-for-data-" + std::to_string(v));
            deleteFile(storagePath + "/proc-ready-" + std::to_string(v));
            deleteFile(storagePath + "/output-status-" + std::to_string(v));
        };
        std::cout << "DONE" << std::endl;


        ////// ALLOCATE SHARED MEMORY FOR CHILD PROCESSES

        // Estimate size of shared memory
        Ciphertext<DCRTPoly> lowLevelCtxt = cryptoContext->Compress(zeroC, 2);
        serializedSize = Serial::serialize(lowLevelCtxt).size();
        sharedMemorySize = serializedSize * numCiphertexts;
        std::cout << "Allocating shared memory (size " << sharedMemorySize
                << " bytes)... " << std::flush;

        // Create the shared memory
        int shm_fd = shm_open("/shared_memory", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
        if (shm_fd == -1)
        {
            perror("shm_open");
            return 1;
        }

        // Resize the shared memory to the required size
        if (ftruncate(shm_fd, sharedMemorySize) == -1)
        {
            perror("ftruncate");
            return 1;
        }

        // Map the shared memory
        sharedMem = static_cast<char*>(mmap(NULL, sharedMemorySize, PROT_READ, MAP_SHARED, shm_fd, 0));
        if (sharedMem == MAP_FAILED)
        {
            perror("mmap");
            return 1;
        }
        std::cout << "DONE" << std::endl;


        ////// PREPARE THE CHILD PROCESSES

        // Spawn subprocesses
        // int numCores = sysconf(_SC_NPROCESSORS_ONLN);  // Number of online processors
        std::cout << "Spawning subprocesses..." << std::endl;
        for (size_t v = 0; v < numCiphertexts; v++)
        {
            pid_t pid = fork();
            if (pid == -1)
            {
                // Fork failed
                std::cerr << "Fork failed" << std::endl;
                return 1;
            }

            if (pid == 0)
            {
                // In child process
                execl(
                    "/usr/bin/numactl",
                    "numactl",
                    ("--cpunodebind=" + std::to_string(v % numNumaNodes)).c_str(),
                    ("--membind=" + std::to_string(v % numNumaNodes)).c_str(),
                    subprocPath.c_str(),
                    std::to_string(v).c_str(),
                    std::to_string(verbose).c_str(),
                    std::to_string(sharedMemorySize).c_str(),
                    std::to_string(serializedSize).c_str(),
                    std::to_string(ndimB).c_str(),
                    std::to_string(indCheby).c_str(),
                    NULL
                );

                // If execl fails
                perror("execl");
                return 1;
            }
            else
            {
                // In parent process
                pidToId[pid] = v; // Store the child's PID -> ciphertext ID
            }
        }

        wait(numCiphertexts, storagePath, "proc-ready-for-data");
        std::cout << "All child processes have spawned." << std::endl;
    }

    double offlineRuntime = 0.0;
    int offlineComm = 0;
    if (partyId != 1)
    {
        offlineRuntime = TOC(t);
        offlineComm = commSize;
        commSize = 0;
    }



    ////////////////////////////////////////////////////////////////////////
    //                     Bob: Datapoints Encryption                     //
    ////////////////////////////////////////////////////////////////////////


    // sync
    if (partyId == 0)
    {
        sendMessage(socket, "");
        receiveMessage(socket);
    }
    else if (partyId == 1)
    {
        receiveMessage(socket);
        sendMessage(socket, "");
    }

    // Bob encrypts his component of the datapoints and sends them to Alice
    std::cout << "Bob encrypts his component of the datapoints xB... " << std::flush;
    if (partyId != 1)
        TIC(t);
    Ciphertext<DCRTPoly> XB[numCompressedCiphertexts][ndimB];
    if (partyId != 0)
    {
        for (size_t d = 0; d < ndimB; d++)
            for (size_t u = 0; u < numCompressedCiphertexts; u++)
            {
                std::vector<double> xBslice(maxNumElCompressed, 0.0);
                size_t i = 0;
                if (k == 2)
                    for (size_t w = 0; w < maxNumEl; w++)
                    {
                        i = u * maxNumElCompressed + w;
                        if (i == n) break;
                        xBslice[w] = xB[i][d];
                    }
                else
                    for (size_t j1 = 0; j1 < k && i < n; j1++)
                        for (size_t j2 = 0; j2 < k && i < n; j2++)
                            for (size_t w = 0; w < maxNumEl; w++)
                            {
                                i = u * maxNumElCompressed + j1 * k * maxNumEl + j2 * maxNumEl + w;
                                if (i == n) break;
                                xBslice[j1 * k * maxNumEl + w * k + j2] = xB[i][d];
                            }
                XB[u][d] = cryptoContext->Encrypt(
                    publicKey,
                    cryptoContext->MakeCKKSPackedPlaintext(xBslice)
                );
            }
    }
    
    if (partyId == -1)
    {
        for (size_t d = 0; d < ndimB; d++)
            for (size_t u = 0; u < numCompressedCiphertexts; u++)
                Serial::serializeToFile(XB[u][d], storagePath + "/XB-"
                                + std::to_string(u) + "-" + std::to_string(d));
    }
    else if (partyId == 1)
    {
        for (size_t d = 0; d < ndimB; d++)
            for (size_t u = 0; u < numCompressedCiphertexts; u++)
                sendMessage(socket, Serial::serialize(XB[u][d]));
    }
    else if (partyId == 0)
    {
        for (size_t d = 0; d < ndimB; d++)
            for (size_t u = 0; u < numCompressedCiphertexts; u++)
            {
                writeToBinFile(
                    storagePath + "/XB-" + std::to_string(u) + "-" + std::to_string(d),
                    receiveMessage(socket));
            }
    }
    std::cout << "DONE" << std::endl;



    ////////////////////////////////////////////////////////////////////////
    //                 Pre-Computations (data-dependent)                  //
    ////////////////////////////////////////////////////////////////////////

    double onlinePrepRuntime = 0.0;
    int onlinePrepComm = 0;
    if (partyId != 1)
    {

        ////// STORING DATAPOINTS FOR CHILD PROCESSES

        std::cout << "Storing datapoints for child processes... " << std::flush;
        Serial::serializeToFile(xA, storagePath + "/xA");
        std::cout << "DONE" << std::endl;


        ////// TRIGGER ALL CHILD PROCESSES TO START PREPARATION

        std::cout << "Child processes: data re-encoding..." << std::flush;
        createFile(storagePath + "/start-preparation");
        wait(numCiphertexts, storagePath, "proc-ready");
        std::cout << "All child processes completed re-encoding." << std::endl;

        onlinePrepRuntime = TOC(t);
        onlinePrepComm = commSize;
        commSize = 0;
    }

    waitForOptimalTemperature(temperatureThreshold);



    ////////////////////////////////////////////////////////////////////////
    //                         LLoyd's Algorithm                          //
    ////////////////////////////////////////////////////////////////////////

    if (partyId != 1)
        TIC(t);

    std::vector<std::vector<std::vector<double>>> centersA;
    std::vector<std::vector<std::vector<double>>> centersB;
    std::vector<double> aliceCompTime;
    for (size_t round = 0; round < numRounds; round++)
    {
        Ciphertext<DCRTPoly> result;
        if (partyId != 1)
        {
            std::cout << std::endl;
            std::cout << "Centers Alice: " << cA << std::endl;
            std::cout << "Centers Bob  : " << cB << std::endl;

            // Alice performs the centers to clusters phase locally
            TIC(t1);
            centersA.push_back(cA);
            centersB.push_back(cB);

            // Storing centers for child processes
            Serial::serializeToFile(cA, storagePath + "/cA");
            Serial::serializeToFile(cB, storagePath + "/cB");

            // Trigger all child processes to start
            createFile(storagePath + "/start-computing");
            sleep(2);
            deleteFile(storagePath + "/start-computing");

            // Generate DP noise (the noise scale is the same for sumA, sumB, and
            // clusterSize)
            std::vector<double> noise = DP::generateGaussianNoise(noiseScale, (k < 3) ? (ndim + 1) * k : 3 * k * maxNumEl);

            // Encrypt the noise as the base variable for aggregating the output
            // from the child processes
            result = cryptoContext->Encrypt(
                publicKey,
                cryptoContext->MakeCKKSPackedPlaintext(noise)
            );

            // Read output of child processes
            std::set<size_t> ongoingProcs; // set of ongoing child processes
            std::vector<std::thread> readOutputThreads;
            for (size_t v = 0; v < numCiphertexts; v++)
                ongoingProcs.insert(v);
            while (ongoingProcs.size() > 0)
            {
                // Track processes that are done computing
                for (auto it = ongoingProcs.begin(); it != ongoingProcs.end(); )
                {
                    size_t v = *it;
                    bool outputIsReady = doesFileExist(storagePath + "/output-status-" + std::to_string(v));
                    if (outputIsReady)
                    {
                        // Remove process from ongoing processes set
                        it = ongoingProcs.erase(it);

                        // Process output in a different thread
                        readOutputThreads.emplace_back(
                            processOutput,
                            std::ref(storagePath),
                            std::ref(result),
                            v,
                            sharedMem,
                            serializedSize
                        );
                    }
                    else it++;
                }
            }
            for (auto& thread : readOutputThreads)
                thread.join(); // Join threads to ensure all of them finish before exiting

            // Store Alice computation time
            aliceCompTime.push_back(TOC(t1));

            // Alice sends the cluster sums and sizes to Bob as a single ciphertext
            size_t encMessageSize = Serial::serialize(result).size();
            std::cout << "Alice sends cluster sums and sizes to Bob (size "
                    << encMessageSize << " bytes)" << std::endl;
        }

        if (partyId == 0)
            sendMessage(socket, Serial::serialize(result));
        else if (partyId == 1)
            result = Serial::deserialize<Ciphertext<DCRTPoly>>(receiveMessage(socket));

        std::string cAserial;
        std::string cBserial;
        if (partyId != 0)
        {
            // Bob decrypts the ciphertext
            cryptoContext->Decrypt(secretKey, result, &plaintext);
            plaintext->SetLength((k < 3) ? (ndim + 1) * k : 3 * k * maxNumEl);
            std::vector<double> resultPtxt = plaintext->GetRealPackedValue();

            // Bob extract sumA, sumB, and clusterSize from the plaintext vector
            std::vector<double> clusterSize(k);
            std::vector<std::vector<double>> sumA(ndimA);
            std::vector<std::vector<double>> sumB(ndimB);
            if (k == 2)
            {
                clusterSize[0] = resultPtxt[0];
                clusterSize[1] = resultPtxt[0];
                size_t index = 1;
                for (size_t d = 0; d < ndimA; d++)
                {
                    sumA[d] = std::vector<double>(
                        resultPtxt.begin() + index,
                        resultPtxt.begin() + index + k
                    );
                    index += k;
                }
                for (size_t d = 0; d < ndimB; d++)
                {
                    sumB[d] = std::vector<double>(
                        resultPtxt.begin() + index,
                        resultPtxt.begin() + index + k
                    );
                    index += k;
                }
            }
            else
            {
                clusterSize.assign(
                    resultPtxt.begin(),
                    resultPtxt.begin() + k
                );
                size_t index = maxNumEl * k;
                for (size_t d = 0; d < ndimA; d++)
                {
                    sumA[d] = std::vector<double>(
                        resultPtxt.begin() + index,
                        resultPtxt.begin() + index + k
                    );
                    index += k;
                }
                index = maxNumEl * 2 * k;
                for (size_t d = 0; d < ndimB; d++)
                {
                    sumB[d] = std::vector<double>(
                        resultPtxt.begin() + index,
                        resultPtxt.begin() + index + k
                    );
                    index += k;
                }
            }

            // Bob computes the new centers
            std::cout << "Bob computes new centers." << std::endl;
            for (size_t j = 0; j < k; j++)
            {
                for (size_t d = 0; d < ndimA; d++)
                    cA[j][d] = sumA[d][j] / clusterSize[j];
                for (size_t d = 0; d < ndimB; d++)
                    cB[j][d] = sumB[d][j] / clusterSize[j];
            }

            // Bob sends the new centers to Alice
            cAserial = Serial::serialize(cA);
            cBserial = Serial::serialize(cB);
            size_t ptxtMessageSize = cAserial.size() + cBserial.size();
            std::cout << "Bob sends updated centers to Alice (size "
                    << ptxtMessageSize << " bytes)." << std::endl;
        }

        if (partyId == 1)
        {
            sendMessage(socket, cAserial);
            sendMessage(socket, cBserial);
        }
        else if (partyId == 0)
        {
            cA = Serial::deserialize<std::vector<std::vector<double>>>(receiveMessage(socket));
            cB = Serial::deserialize<std::vector<std::vector<double>>>(receiveMessage(socket));
        }

        if (partyId != 1)
        {
            double ctxtLoss = computeLoss(xA, xB, cA, cB);
            double ptxtLoss = computeLoss(xA, xB, expectedCenters[0], expectedCenters[1]);
            std::cout << "K-means loss of encrypted approach: " << ctxtLoss << std::endl;
            std::cout << "K-means loss of plaintext approach: " << ptxtLoss << std::endl;
            std::cout << std::endl;
        }

        // Checking for convergence
        if (partyId != 1)
            if (numRounds == -1UL)
            {
                std::vector<size_t> clusterSize(k);
                for (size_t i = 0; i < n; i++)
                {
                    size_t minIndex = 0;
                    double minDist = ndim;
                    for (size_t j = 0; j < k; j++)
                    {
                        double currDist = 0.0;
                        for (size_t d = 0; d < ndimA; d++)
                            currDist += std::pow(cA[j][d] - xA[i][d], 2);
                        for (size_t d = 0; d < ndimB; d++)
                            currDist += std::pow(cB[j][d] - xB[i][d], 2);
                        if (currDist < minDist)
                        {
                            minIndex = j;
                            minDist = currDist;
                        }
                    }
                    clusterSize[minIndex] += 1;
                }
                bool hasConverged = true;
                for (size_t it = MAX(0, (int) centersA.size() - (int) convergenceWindowSize); it < centersA.size() && hasConverged; it++)
                    for (size_t j = 0; j < k && hasConverged; j++)
                    {
                        double dist = 0;
                        for (size_t d = 0; d < ndimA; d++)
                            dist += std::pow(centersA[it][j][d] - cA[j][d], 2);
                        for (size_t d = 0; d < ndimB; d++)
                            dist += std::pow(centersB[it][j][d] - cB[j][d], 2);
                        if (dist > std::pow(convergenceRadius / clusterSize[j], 2))
                            hasConverged = false;
                    }
                if (hasConverged)
                    break;
            }
    }

    double onlineRuntime;
    int onlineComm;
    if (partyId != 1)
    {
        onlineRuntime = TOC(t);
        onlineComm = commSize;
    }

    // Print result
    if (partyId != 1)
    {
        std::cout << std::endl;
        std::cout << "Computed cA: " << cA << std::endl;
        std::cout << "Computed cB: " << cB << std::endl;
        std::cout << "Expected cA: " << expectedCenters[0] << std::endl;
        std::cout << "Expected cB: " << expectedCenters[1] << std::endl;
        std::cout << std::endl;
    }



    ////////////////////////////////////////////////////////////////////////
    //                    Accuracy and Runtime Metrics                    //
    ////////////////////////////////////////////////////////////////////////

    if (partyId != 1)
    {
        // Compute l2 distance between computed and expected centroids
        double maxDist = 0.0;
        double avgDist = 0.0;
        std::vector<double> cerror(k, 0.0);
        for (size_t j = 0; j < k; j++)
        {
            for (size_t d = 0; d < ndimA; d++)
                cerror[j] += std::pow(cA[j][d] - expectedCenters[0][j][d], 2);
            for (size_t d = 0; d < ndimB; d++)
                cerror[j] += std::pow(cB[j][d] - expectedCenters[1][j][d], 2);
            cerror[j] = std::sqrt(cerror[j]);
            maxDist = MAX(maxDist, cerror[j]);
            avgDist += cerror[j];
        }
        avgDist /= k;
        std::cout << "Centroid error: " << cerror << std::endl;
        std::cout << "           max: " << maxDist << std::endl;
        std::cout << "           avg: " << avgDist << std::endl;
        std::cout << std::endl;

        // Compute k-means loss
        double ctxtLoss = computeLoss(xA, xB, cA, cB);
        double ptxtLoss = computeLoss(xA, xB, expectedCenters[0], expectedCenters[1]);
        std::cout << "K-means loss of encrypted approach: " << ctxtLoss << std::endl;
        std::cout << "K-means loss of plaintext approach: " << ptxtLoss << std::endl;
        std::cout << std::endl;

        // Compute clusters accuracy
        double ctxtAcc = -1.0;
        double ptxtAcc = -1.0;
        if (testClusters)
        {
            std::cout << "Cluster accuracy of encrypted approach" << std::endl;
            ctxtAcc = clusterAccuracy(xA, xB, cA, cB, groundTruth);
            std::cout << "Cluster accuracy of plaintext approach" << std::endl;
            ptxtAcc = clusterAccuracy(xA, xB, expectedCenters[0], expectedCenters[1], groundTruth);
        }

        // Compute average computation time
        double avgtime = 0.0;
        for (double t : aliceCompTime)
            avgtime += t;
        avgtime /= aliceCompTime.size();

        // Print runtime and communication
        int comm = offlineComm + onlinePrepComm + onlineComm;
        double runtime = offlineRuntime + onlinePrepRuntime + onlineRuntime;
        std::cout << "                 comm. (bytes)     runtime (s)" << std::endl;
        std::cout << "Offline phase:" << std::setw(16) << offlineComm       << std::setw(16) << offlineRuntime      << std::endl;
        std::cout << " Online prep.:" << std::setw(16) << onlinePrepComm    << std::setw(16) << onlinePrepRuntime   << std::endl;
        std::cout << " Online Lloyd:" << std::setw(16) << onlineComm        << std::setw(16) << onlineRuntime       << std::endl;
        std::cout << "          TOT:" << std::setw(16) << comm              << std::setw(16) << runtime             << std::endl;
        std::cout << std::endl;
        std::cout << "Average computation time for one round: " << avgtime << "s" << std::endl;
        std::cout << std::endl;

        // Store accuracy, communication, and runtime metrics to file
        std::ofstream file;
        file.open("metrics.csv", std::ios_base::app);  // Open in append mode
        if (file.is_open())
        {
            file           << n
                    << "," << k
                    << "," << ndimA
                    << "," << ndimB
                    << "," << numRounds
                    << "," << epsilonTotal
                    << "," << deltaTotal
                    << "," << numCiphertexts
                    << "," << maxDist
                    << "," << avgDist
                    << "," << ctxtLoss
                    << "," << ptxtLoss
                    << "," << ctxtAcc
                    << "," << ptxtAcc
                    << "," << offlineComm
                    << "," << onlinePrepComm
                    << "," << onlineComm
                    << "," << offlineRuntime
                    << "," << onlinePrepRuntime
                    << "," << onlineRuntime
                    << "," << avgtime
                    << std::endl;
            file.close();  // Close the file
        }
        else
            std::cerr << "Unable to open file: metrics.csv" << std::endl;
    }



    ////////////////////////////////////////////////////////////////////////
    //                        Closing and Cleaning                        //
    ////////////////////////////////////////////////////////////////////////

    // Stop all child processes
    for (size_t v = 0; v < numCiphertexts; v++)
        createFile(storagePath + "/stop");

    // Wait for all child processes to finish
    for (auto& pid : pidToId) {
        int status;
        waitpid(pid.first, &status, 0);
    }

    // Detach and close shared memory
    munmap(sharedMem, sharedMemorySize);
    shm_unlink("/shared_memory");

    // Close sockets
    if (partyId == 0)
        close(server_fd);
    if (partyId == 0 || partyId == 1)
        close(socket);

    return 0;

}
