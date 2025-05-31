#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-ptxt.h"
#include "utils-matrices.h"

#include "math/chebyshev.h"

#include <filesystem>


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


/**
 * @brief Change-point detection using triplet comparison.
 * 
 * This function detects change points in a time series by comparing triplets
 * of consecutive differences. It divides the time series into blocks of a
 * specified size (default is the square root of the time series length),
 * counts how many triplets in each block are not monotonically increasing or
 * decreasing, and computes a CUSUM statistic to find the change point.
 * 
 * @param x The input time series as a vector of doubles.
 * @param blockSize The size of the blocks to divide the time series into
 *        (default is 0).
 * @return The index of the detected change point in the time series.
 * 
 * @note If `blockSize` is set to 0, it defaults to the square root of the
 * length of the time series.
 */
size_t cpd
(
    const std::vector<double> x,
    size_t blockSize = 0,
    const int changePointType = 0  // 0 - mean, 1 - variance, 2 - distribution
)
{
    // Set default block size to sqrt(x.size()) if not provided
    size_t n = x.size();
    if (blockSize == 0)
        blockSize = static_cast<size_t>(std::sqrt(n));
    const size_t numBlocks = n / blockSize;

    size_t changePoint = -1;

    if (changePointType == 0)
    {
        std::cout << "Change-point type: mean" << std::endl;

        // For each block, sum all its elements.
        std::vector<double> q(numBlocks, 0.0);
        for (size_t i = 0; i < numBlocks; i++)
        {
            // Initialize the count for the current block
            q[i] = 0;

            // Determine the start and end indices for the current block
            size_t start = i * blockSize;
            size_t end = std::min(start + blockSize, n - 1);

            // Sum the elements in the current block
            for (size_t j = start; j < end; j++)
                q[i] += x[j];
        }

        // CUSUM statistic: for each i = 1, ..., numBlocks - 1
        // compute s[i] = | sum(q[0:i]) - i/numBlocks * sum(q[0:numBlocks]) |
        std::vector<double> s(numBlocks - 1, 0.0);
        double sumQ = std::accumulate(q.begin(), q.end(), 0.0);
        for (size_t i = 0; i < numBlocks - 1; i++)
        {
            double sumQ_i = std::accumulate(q.begin(), q.begin() + i + 1, 0.0);
            s[i] = std::abs(sumQ_i - (i + 1) * sumQ / numBlocks);
        }

        // Find the index of the maximum value in s
        size_t maxIndex = 0;
        double maxValue = s[0];
        for (size_t i = 1; i < s.size(); i++)
        {
            if (s[i] > maxValue)
            {
                maxValue = s[i];
                maxIndex = i;
            }
        }

        // Adjust the index to be relative to the original vector x
        changePoint = (maxIndex + 1) * blockSize;
        if (changePoint >= n)
            changePoint = n - 1; // Ensure the change point is within bounds
    }
    else if (changePointType == 1)
    {
        std::cout << "Change-point type: variance" << std::endl;

        // For each block, compute its variance.
        std::vector<double> q(numBlocks, 0.0);
        for (size_t i = 0; i < numBlocks; i++)
        {
            // Initialize the count for the current block
            q[i] = 0;
            double mean = 0.0;

            // Determine the start and end indices for the current block
            size_t start = i * blockSize;
            size_t end = std::min(start + blockSize, n - 1);

            // Compute the mean of the current block
            for (size_t j = start; j < end; j++)
                mean += x[j];
            mean /= (end - start);

            // Compute the variance of the current block
            for (size_t j = start; j < end; j++)
                q[i] += (x[j] - mean) * (x[j] - mean);
        }

        // CUSUM statistic: for each i = 1, ..., numBlocks - 1
        // compute s[i] = | sum(q[0:i]) - i/numBlocks * sum(q[0:numBlocks]) |
        std::vector<double> s(numBlocks - 1, 0.0);
        double sumQ = std::accumulate(q.begin(), q.end(), 0.0);
        for (size_t i = 0; i < numBlocks - 1; i++)
        {
            double sumQ_i = std::accumulate(q.begin(), q.begin() + i + 1, 0.0);
            s[i] = std::abs(sumQ_i - (i + 1) * sumQ / numBlocks);
        }

        // Find the index of the maximum value in s
        size_t maxIndex = 0;
        double maxValue = s[0];
        for (size_t i = 1; i < s.size(); i++)
        {
            if (s[i] > maxValue)
            {
                maxValue = s[i];
                maxIndex = i;
            }
        }

        // Adjust the index to be relative to the original vector x
        changePoint = (maxIndex + 1) * blockSize;
        if (changePoint >= n)
            changePoint = n - 1; // Ensure the change point is within bounds
    }
    else if (changePointType == 2)
    {
        std::cout << "Change-point type: distribution" << std::endl;

        // For each block, iterate over the subsequent triplets and count how many
        // of them are not monotonically increasing or decreasing.
        std::vector<double> q(numBlocks, 0.0);
        for (size_t i = 0; i < numBlocks; i++)
        {
            // Initialize the count for the current block
            q[i] = 0;

            // Determine the start and end indices for the current block
            size_t start = i * blockSize;
            size_t end = std::min(start + blockSize, n - 1);

            // Iterate over the triplets in the current block
            for (size_t j = start; j < end - 2; j++)
            {
                // Check if the triplet is not monotonically increasing or decreasing
                if (!((x[j] < x[j + 1] && x[j + 1] < x[j + 2]) ||
                    (x[j] > x[j + 1] && x[j + 1] > x[j + 2])))
                {
                    q[i]++;
                }
            }

            // Normalize by the number of triplets
            std::cout << "q[" << i << "] = " << (size_t) q[i] << std::endl;
            q[i] /= blockSize - 2;
        }

        // CUSUM statistic: for each i = 1, ..., numBlocks - 1
        // compute s[i] = | sum(q[0:i]) - i/numBlocks * sum(q[0:numBlocks]) |
        std::vector<double> s(numBlocks - 1, 0.0);
        double sumQ = std::accumulate(q.begin(), q.end(), 0.0);
        for (size_t i = 0; i < numBlocks - 1; i++)
        {
            double sumQ_i = std::accumulate(q.begin(), q.begin() + i + 1, 0.0);
            s[i] = std::abs(sumQ_i - (i + 1) * sumQ / numBlocks);
        }

        // Find the index of the maximum value in s
        size_t maxIndex = 0;
        double maxValue = s[0];
        for (size_t i = 1; i < s.size(); i++)
        {
            if (s[i] > maxValue)
            {
                maxValue = s[i];
                maxIndex = i;
            }
        }

        // Adjust the index to be relative to the original vector x
        changePoint = (maxIndex + 1) * blockSize;
        if (changePoint >= n)
            changePoint = n - 1; // Ensure the change point is within bounds
    }
    else
    {
        throw std::runtime_error("Unknown change-point type");
    }
    
    return changePoint;
}


/**
 * @brief Load a dataset from a CSV file.
 * 
 * This function reads a single line from the specified CSV file, parses the
 * comma-separated values into a vector of doubles, and optionally extracts a
 * label if provided. The data is normalized to the range [0, 1] if specified.
 * 
 * @param filepath Path to the CSV file containing the dataset.
 * @param[out] x Vector to store the parsed time series data.
 * @param[out] label Pointer to an integer to store the label,if provided
 *             (default is nullptr).
 * @param normalize Boolean flag to indicate whether to normalize the data
 *        (default is true).
 * @param verbose Boolean flag to indicate whether to print the loaded data
 *        (default is false).
 * 
 * @throws std::runtime_error If the file cannot be opened, or if there are
 * invalid numbers in the file, or if the label is invalid.
 * 
 * @note The function assumes that the CSV file contains a single line with
 * comma-separated values, where the last value is the label if `label` is not
 * null.
 */
void loadDataset
(
    const std::string &filepath,
    std::vector<double> &x,
    int *label=nullptr,
    const bool normalize=true,
    const bool verbose=false
)
{
    // read time series from .csv file
    // there should be only one line with comma-separated values, one for each time point
    // if label is provided, it will be read as the last value in the line

    // Opening the file
    std::ifstream infile(filepath);
    if (!infile)
        throw std::runtime_error("Error: Unable to open file " + filepath);
    
    // Reading the first line
    std::string line;
    if (!std::getline(infile, line))
    {
        infile.close();
        throw std::runtime_error("Error: Unable to read from file " + filepath);
    }
    infile.close();

    // Splitting the line by comma
    std::stringstream ss(line);
    std::string value;
    while (std::getline(ss, value, ','))
    {
        try
        {
            x.push_back(std::stod(value));
        }
        catch (const std::invalid_argument &e)
        {
            throw std::runtime_error("Error: Invalid number in file " + filepath);
        }
    }

    // If label is provided, read it as the last value in the line
    if (label != nullptr)
    {
        if (x.empty())
            throw std::runtime_error("Error: No data points found in file " + filepath);
        
        // Assuming the label is the last value in the line
        try
        {
            *label = static_cast<int>(x.back());
            x.pop_back(); // Remove label from data points
        }
        catch (const std::invalid_argument &e)
        {
            throw std::runtime_error("Error: Invalid label in file " + filepath);
        }
    }

    // Check if we have at least one data point
    if (x.empty())
        throw std::runtime_error("Error: No data points found in file " + filepath);
    
    // Check if the label is valid
    if (label != nullptr && (*label < 0 || *label >= static_cast<int>(x.size())))
        throw std::runtime_error("Error: Invalid label in file " + filepath);

    // Print the loaded data (first 10 values and label) if verbose is true
    if (verbose)
    {
        std::cout << "Time series with " << x.size() << " points loaded from " 
                  << filepath << ":" << std::endl;
        for (size_t i = 0; i < std::min(x.size(), size_t(10)); ++i)
            std::cout << x[i] << " ";
        if (x.size() > 10)
            std::cout << "...   ";
        if (label != nullptr)
            std::cout << "[Label: " << *label << "]" << std::endl;
        else
            std::cout << "[No label provided]" << std::endl;
    }

    // Normalize the data if required
    if (normalize)
    {
        double minVal = *std::min_element(x.begin(), x.end());
        double maxVal = *std::max_element(x.begin(), x.end());
        if (maxVal - minVal > 0.0)
        {
            for (double &val : x)
                val = (val - minVal) / (maxVal - minVal);
        }
        else
        {
            // If all values are the same, normalize to 0.5
            std::fill(x.begin(), x.end(), 0.5);
        }
    }

}


/**
 * @brief Print the contents of a ciphertext.
 * 
 * This function decrypts a ciphertext using the provided secret key and prints
 * the decrypted values in a formatted manner.
 * 
 * @param cryptoContext The CryptoContext used for decryption.
 * @param secretKey The secret key used for decryption.
 * @param ctxt The ciphertext to be decrypted and printed.
 * @param label An optional label to be printed before the ciphertext values
 *        (default is an empty string).
 * @param numCols The number of columns to print in the output (default is 10).
 * @param numRows The number of rows to print in the output (default is 1).
 */
void printCtxt
(
    const CryptoContext<DCRTPoly> &cryptoContext,
    const PrivateKey<DCRTPoly> &secretKey,
    const Ciphertext<DCRTPoly> &ctxt,
    const std::string &label = "",
    const size_t numCols=10,
    const size_t numRows=1,
    const size_t maxSize=20,
    const bool dumpToFile=false
)
{
    Plaintext plaintext;
    cryptoContext->Decrypt(secretKey, ctxt, &plaintext);
    std::vector<double> data = plaintext->GetRealPackedValue();

    std::cout << std::endl
              << label << " (level: " << ctxt->GetLevel() << ") [min-max: "
              << *std::min_element(data.begin(), data.end()) << " - "
              << *std::max_element(data.begin(), data.end()) << "]:"
              << std::endl;
    for (size_t i = 0; i < MIN(numRows, maxSize); ++i)
    {
        for (size_t j = 0; j < MIN(numCols, maxSize) && i * numCols + j < data.size(); ++j)
        {
            std::cout << std::fixed << std::setprecision(3) 
                      << data[i * numCols + j] << " ";
        }
        std::cout << std::endl;
    }

    // store matrix as csv in a file data/out/<label>.csv
    if (dumpToFile && !label.empty())
    {
        // Ensure the output directory exists
        std::filesystem::create_directories("data/out");

        std::ofstream outFile("data/out/" + label + ".csv");
        if (!outFile)
        {
            std::cerr << "Error: Unable to open output file for writing." << std::endl;
            return;
        }
        for (size_t i = 0; i < numRows; ++i)
        {
            for (size_t j = 0; j < numCols && i * numCols + j < data.size(); ++j)
            {
                outFile << std::fixed << std::setprecision(3) 
                        << data[i * numCols + j];
                if (j < numCols - 1)
                    outFile << ",";
            }
            outFile << std::endl;
        }
        outFile.close();
    }
}



int main(int argc, char* argv[])
{

    std::cout << "Change-point detection on encrypted time series"
        << std::endl << std::endl;

    bool verbose = false;

    // File with datapoints
    std::string filepath = "";

    // Block size
    size_t blockSize = 0;

    // Approximation degrees
    usint compareTDepth = 10;
    usint compareADepth = 10;

    // Is the label provided
    bool isLabeled = false;

    // Change-point type
    // 0 - mean, 1 - variance, 2 - distribution
    int changePointType = 0;

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
        else if (arg == "-b" || arg == "--block-size")
            if (i + 1 < argc)
                blockSize = std::stoul(argv[++i]);
            else
            {
                std::cerr << "Error: --block-size requires a size argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-ctd" || arg == "--cmpt-degree")
            if (i + 1 < argc)
                compareTDepth = std::stoi(argv[++i]);
            else
            {
                std::cerr << "Error: --cmp-degree requires a degree argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-cad" || arg == "--cmpa-degree")
            if (i + 1 < argc)
                compareADepth = std::stoi(argv[++i]);
            else
            {
                std::cerr << "Error: --cmp-degree requires a degree argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-l" || arg == "--labeled")
            isLabeled = true;
        else if (arg == "-t" || arg == "--type")
            if (i + 1 < argc)
            {
                std::string type = argv[++i];
                if (type == "mean")
                    changePointType = 0;
                else if (type == "variance")
                    changePointType = 1;
                else if (type == "distribution")
                    changePointType = 2;
                else
                {
                    std::cerr << "Error: Unknown change-point type: " << type
                              << ". Supported types are: mean, variance, distribution."
                              << std::endl;
                    return 1;
                }
            }
            else
            {
                std::cerr << "Error: --type requires a type argument."
                          << std::endl;
                return 1;
            }
        else if (arg == "-h" || arg == "--help")
        {
            std::cout << "Usage: clustering "
                      << "-f|--file <filepath> "
                      << "[-ctd|--cmpt-degree <degree>] "
                      << "[-cad|--cmpa-degree <degree>] "
                      << "[-l|--labeled] "
                      << "[-v|--verbose] "
                      << "[-h|--help] "
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
    if (filepath.length() == 0)
    {
        std::cerr << "Missing datapoints filepath." << std::endl;
        return 1;
    }

    std::vector<double> x;  // time series data
    int label;              // label of the time series (if provided)

    loadDataset(
        filepath,
        x,
        isLabeled ? &label : nullptr,
        true,  // normalize
        verbose
    );

    // Setting default block size if not provided
    blockSize = (blockSize == 0) ? static_cast<size_t>(std::sqrt(x.size())) : blockSize;

    const size_t n = x.size();  // number of time points
    if (blockSize == 0)
        blockSize = static_cast<size_t>(std::sqrt(n));
    const size_t numBlocks = n / blockSize;
    std::cout << "Change-point type: " 
              << (changePointType == 0 ? "mean" : changePointType == 1 ? "variance" : "distribution")
              << std::endl;
    std::cout << "Number of time points: " << n << std::endl;
    std::cout << "Block size: " << blockSize << std::endl;
    std::cout << "Number of blocks: " << numBlocks << std::endl;
    std::cout << "Triplets comparison depth: " << compareTDepth << std::endl;
    std::cout << "Argmax comparison depth: " << compareADepth << std::endl;
    std::cout << std::endl;

    // num columns is the next power of 2 greater than or equal to blockSize
    const size_t numColumns =
        static_cast<size_t>(std::pow(2, std::ceil(std::log2(blockSize))));

    // Initialize timers
    auto t = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSeconds;
    std::cout << std::fixed << std::setprecision(3);
    #undef TOC
    #define TOC(t) (elapsedSeconds = timeNow() - t).count()



    ////////////////////////////////////////////////////////////////////////
    //                          CPD in Plaintext                          //
    ////////////////////////////////////////////////////////////////////////

    size_t expectedChangePoint = cpd(x, blockSize, changePointType);
    std::cout << "Expected change point: " << expectedChangePoint << std::endl;



    ////////////////////////////////////////////////////////////////////////
    //                       Setting up CKKS scheme                       //
    ////////////////////////////////////////////////////////////////////////

    bool cmpAdv               = true;
    usint dg_ct               = 4;
    usint df_ct               = 2;
    usint dg_ca               = 4;
    usint df_ca               = 2;
    if (changePointType < 2)
    {
        compareTDepth = 0;
        dg_ct = 0;
        df_ct = 0;
    }
    const usint integralPrecision   = std::floor(std::log2(n)) + 1;
    const usint decimalPrecision    = 60 - (int) integralPrecision;
    const usint cmpCost             = cmpAdv ? 4 * (dg_ct + df_ct + dg_ca + df_ca) + 2 : compareTDepth + compareADepth;
    const usint multiplicativeDepth = cmpCost + 6 + LOG2(numColumns)
                                      + (changePointType == 1 ? 3 : 0)
                                      + (changePointType == 2 ? 2 : 0);
    const usint numSlots            = numColumns * numColumns;

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
        integralPrecision,
        decimalPrecision,
        multiplicativeDepth,
        numSlots,
        false, 0, true
    );

    // Generating public/private key pair, relinearization, and rotation keys
    // << positive
    // >> negative
    // assuming k <= n
    std::vector<int32_t> indices;
    for (size_t i = 0; i < LOG2(numColumns); i++)
    {
        indices.push_back(1 << i);
        indices.push_back(-(1 << i));
        // indices.push_back(1 << (LOG2(numColumns) + i));
        indices.push_back(-(1 << (LOG2(numColumns) + i)));
        indices.push_back(-(numColumns * (numColumns - 1) / (1 << (i + 1))));
    }
    
    KeyPair<DCRTPoly> keyPair = keyGeneration(
        cryptoContext,
        indices,
        numSlots,
        false, true
    );



    ////////////////////////////////////////////////////////////////////////
    //                       Time Series Encryption                       //
    ////////////////////////////////////////////////////////////////////////


    // Padding the time series data to fit the number of slots. This can be
    // avoided by using an additional O(log(N)) rotations during the algorithm.
    std::vector<double> paddedX(numSlots, 0.0);
    for (size_t i = 0; i < numBlocks; i++)
        for (size_t j = 0; j < blockSize; j++)
            paddedX[i * numColumns + j] = x[i * blockSize + j];

    // Encrypting the time series data
    Ciphertext<DCRTPoly> X = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(paddedX)
    );
    printCtxt(cryptoContext, keyPair.secretKey, X, "Encrypted time series", numColumns, numBlocks);



    ////////////////////////////////////////////////////////////////////////
    //                Pre-Computations (data-independent)                 //
    ////////////////////////////////////////////////////////////////////////


    std::cout << "Generating coefficients for triplet comparison...";

    std::vector<double> cmpTCoeffs = generateCompareCoeffs(
        -1.0, 1.0, depth2degree(compareTDepth));
    std::vector<double> cmpACoeffs = generateCompareCoeffs(
        -1.0, 1.0, depth2degree(compareADepth));
    std::vector<double> rectMask(numSlots, 0.0);
    for (size_t i = 0; i < numBlocks; i++)
        for (size_t j = 0; j < blockSize - 2; j++)
            rectMask[i * numColumns + j] = 1.0;
    std::vector<double> triangularMask(numSlots, 0.0);
    for (size_t i = 0; i < numBlocks; i++)
        for (size_t j = i; j < numBlocks; j++)
            triangularMask[i * numColumns + j] = 1.0;
    std::vector<double> multConstant(numSlots, 0.0);
    for (size_t i = 0; i < numBlocks; i++)
        multConstant[i] = 1.0 * (i + 1) / numBlocks;
    std::vector<double> diagonal(numSlots, 0.0);
    for (size_t i = 0; i < numBlocks; i++)
        diagonal[i * numColumns + i] = 0.5;
    for (size_t i = numBlocks; i < numColumns; i++)
        for (size_t j = 0; j < numBlocks; j++)
            diagonal[i * numColumns + j] = 1.0;
    std::vector<double> rectMask2(numSlots, 0.0);
    for (size_t i = 0; i < numBlocks; i++)
        for (size_t j = 0; j < numBlocks; j++)
            rectMask2[i * numColumns + j] = 1.0;

    std::cout << "DONE" << std::endl;



    ////////////////////////////////////////////////////////////////////////
    //                                CPD                                 //
    ////////////////////////////////////////////////////////////////////////


    std::cout << "Starting change-point detection..." << std::endl;
    TIC(t);

    Ciphertext<DCRTPoly> R;

    if (changePointType == 0)
    {
        std::cout << "Change-point type: mean" << std::endl;

        Ciphertext<DCRTPoly> Q = sumColumns(X, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Sum of blocks", numColumns, numBlocks);

        Q = (1.0 / blockSize) * Q;
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Mean of blocks", numColumns, numBlocks);

        Ciphertext<DCRTPoly> QT = replicateColumn(sumRows(Q, numColumns, true), numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, QT, "QT", numColumns, numBlocks);

        Ciphertext<DCRTPoly> QP = replicateColumn(Q, numColumns) * triangularMask;
        QP = sumRows(QP, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, QP, "QP", numColumns, numBlocks);

        Ciphertext<DCRTPoly> S = QP - multConstant * QT;
        S = S * S;
        printCtxt(cryptoContext, keyPair.secretKey, S, "CUSUM statistic", numColumns, numBlocks);

        // Compute argmax of S
        S = S * (2.0 / numBlocks);
        Ciphertext<DCRTPoly> SR = replicateRow(S, numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, SR, "SR", numColumns, numBlocks);

        Ciphertext<DCRTPoly> SC = replicateColumn(transposeRow(S, numColumns, true), numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, SC, "SC", numColumns, numBlocks);

        Ciphertext<DCRTPoly> CMP = cmpAdv ?
            compareAdv(SR, SC, dg_ca, df_ca) :
            cryptoContext->EvalChebyshevSeries(SR - SC, cmpACoeffs, -1.0, 1.0);
        CMP = CMP * rectMask2 + diagonal;
        printCtxt(cryptoContext, keyPair.secretKey, CMP, "Argmax comparison result", numColumns, numColumns);

        // Multiply the rows of CMP together
        for (size_t i = 0; i < LOG2(numColumns); i++)
        {
            CMP = CMP * (CMP >> (1 << (LOG2(numColumns) + i)));
            printCtxt(cryptoContext, keyPair.secretKey, CMP, "Argmax comparison result after multiplication " + std::to_string(i + 1), numColumns, numColumns);
        }
        R = CMP;    // argmax of S
        printCtxt(cryptoContext, keyPair.secretKey, R, "Computed change-point", numColumns, numBlocks);
    }
    else if (changePointType == 1)
    {
        std::cout << "Change-point type: variance" << std::endl;

        Ciphertext<DCRTPoly> Q = sumColumns(X, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Sum of blocks", numColumns, numBlocks);

        Q = (1.0 / blockSize) * Q;
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Mean of blocks", numColumns, numBlocks);

        Q = replicateColumn(Q, numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Replicated mean of blocks", numColumns, numBlocks);

        Q = X - Q;
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Difference from mean", numColumns, numBlocks);

        Q = (1.0 / (blockSize - 1)) * sumColumns(Q * Q, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Variance of blocks", numColumns, numBlocks);

        Ciphertext<DCRTPoly> QT = replicateColumn(sumRows(Q, numColumns, true), numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, QT, "QT", numColumns, numBlocks);

        Ciphertext<DCRTPoly> QP = replicateColumn(Q, numColumns) * triangularMask;
        QP = sumRows(QP, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, QP, "QP", numColumns, numBlocks);

        Ciphertext<DCRTPoly> S = QP - multConstant * QT;
        S = S * S;
        printCtxt(cryptoContext, keyPair.secretKey, S, "CUSUM statistic", numColumns, numBlocks);

        // Compute argmax of S
        S = S * (2.0 / numBlocks);
        Ciphertext<DCRTPoly> SR = replicateRow(S, numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, SR, "SR", numColumns, numBlocks);

        Ciphertext<DCRTPoly> SC = replicateColumn(transposeRow(S, numColumns, true), numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, SC, "SC", numColumns, numBlocks);

        Ciphertext<DCRTPoly> CMP = cmpAdv ?
            compareAdv(SR, SC, dg_ca, df_ca) :
            cryptoContext->EvalChebyshevSeries(SR - SC, cmpACoeffs, -1.0, 1.0);
        CMP = CMP * rectMask2 + diagonal;
        printCtxt(cryptoContext, keyPair.secretKey, CMP, "Argmax comparison result", numColumns, numColumns);

        // Multiply the rows of CMP together
        for (size_t i = 0; i < LOG2(numColumns); i++)
        {
            CMP = CMP * (CMP >> (1 << (LOG2(numColumns) + i)));
            printCtxt(cryptoContext, keyPair.secretKey, CMP, "Argmax comparison result after multiplication " + std::to_string(i + 1), numColumns, numColumns);
        }
        R = CMP;    // argmax of S
        printCtxt(cryptoContext, keyPair.secretKey, R, "Computed change-point", numColumns, numBlocks);
    }
    else if (changePointType == 2)
    {
        std::cout << "Change-point type: distribution" << std::endl;
        
        Ciphertext<DCRTPoly> C = cmpAdv ?
            compareAdv(X, (X << 1), dg_ct, df_ct) :
            cryptoContext->EvalChebyshevSeries(X - (X << 1), cmpTCoeffs, -1.0, 1.0);

        printCtxt(cryptoContext, keyPair.secretKey, C, "Triplet comparison result", numColumns, numBlocks);

        Ciphertext<DCRTPoly> D = C + (C << 1);
        printCtxt(cryptoContext, keyPair.secretKey, D, "Triplet comparison sum", numColumns, numBlocks);

        Ciphertext<DCRTPoly> E = D - 1;
        Ciphertext<DCRTPoly> F = 1 - (E * E);
        printCtxt(cryptoContext, keyPair.secretKey, E, "Triplet comparison sum - 1", numColumns, numBlocks);
        printCtxt(cryptoContext, keyPair.secretKey, F, "Triplet comparison sum squared", numColumns, numBlocks);

        Ciphertext<DCRTPoly> Q = sumColumns(F * rectMask, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Good triplets counter", numColumns, numBlocks);

        Q = (1.0 / (blockSize - 2)) * Q;
        printCtxt(cryptoContext, keyPair.secretKey, Q, "Good triplets counter normalized", numColumns, numBlocks);

        Ciphertext<DCRTPoly> QT = replicateColumn(sumRows(Q, numColumns, true), numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, QT, "QT", numColumns, numBlocks);

        Ciphertext<DCRTPoly> QP = replicateColumn(Q, numColumns) * triangularMask;
        QP = sumRows(QP, numColumns, true);
        printCtxt(cryptoContext, keyPair.secretKey, QP, "QP", numColumns, numBlocks);

        Ciphertext<DCRTPoly> S = QP - multConstant * QT;
        S = S * S;
        printCtxt(cryptoContext, keyPair.secretKey, S, "CUSUM statistic", numColumns, numBlocks);

        // Compute argmax of S
        S = S * (2.0 / numBlocks);
        Ciphertext<DCRTPoly> SR = replicateRow(S, numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, SR, "SR", numColumns, numBlocks);

        Ciphertext<DCRTPoly> SC = replicateColumn(transposeRow(S, numColumns, true), numColumns);
        printCtxt(cryptoContext, keyPair.secretKey, SC, "SC", numColumns, numBlocks);

        Ciphertext<DCRTPoly> CMP = cmpAdv ?
            compareAdv(SR, SC, dg_ca, df_ca) :
            cryptoContext->EvalChebyshevSeries(SR - SC, cmpACoeffs, -1.0, 1.0);
        CMP = CMP * rectMask2 + diagonal;
        printCtxt(cryptoContext, keyPair.secretKey, CMP, "Argmax comparison result", numColumns, numColumns);

        // Multiply the rows of CMP together
        for (size_t i = 0; i < LOG2(numColumns); i++)
        {
            CMP = CMP * (CMP >> (1 << (LOG2(numColumns) + i)));
            printCtxt(cryptoContext, keyPair.secretKey, CMP, "Argmax comparison result after multiplication " + std::to_string(i + 1), numColumns, numColumns);
        }
        R = CMP;    // argmax of S
        printCtxt(cryptoContext, keyPair.secretKey, R, "Computed change-point", numColumns, numBlocks);
    }
    else
    {
        std::cerr << "Unknown change-point type: " << changePointType << std::endl;
        return 1;
    }

    double runtime = TOC(t);
    std::cout << "CPD runtime: " << runtime << "s" << std::endl;



    ////////////////////////////////////////////////////////////////////////
    //                          Printing Result                           //
    ////////////////////////////////////////////////////////////////////////


    Plaintext plaintext;
    cryptoContext->Decrypt(keyPair.secretKey, R, &plaintext);
    plaintext->SetLength(blockSize);

    // Extract the change-point index from the plaintext as argmax of the vector
    std::vector<double> argmax = plaintext->GetRealPackedValue();
    size_t maxIndex = 0;
    double maxVal = argmax[0];
    for (size_t i = 1; i < argmax.size(); i++)
    {
        if (argmax[i] > maxVal)
        {
            maxVal = argmax[i];
            maxIndex = i;
        }
    }
    size_t changePoint = (maxIndex + 1) * blockSize;
    std::cout << "Computed change-point: " << changePoint << std::endl;
    std::cout << "Expected change-point: " << expectedChangePoint << std::endl;


    return 0;

}
