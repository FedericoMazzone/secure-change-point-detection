#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-ptxt.h"
#include "serial.h"

#include "math/chebyshev.h"

#include <cassert>
#include <fcntl.h>      // For O_RDWR
#include <filesystem>
#include <sys/mman.h>   // For shm_open, mmap, munmap, MAP_SHARED, PROT_READ, PROT_WRITE, MAP_FAILED
#include <sys/stat.h>   // For S_IRUSR, S_IWUSR
#include <unistd.h>     // For close


using namespace lbcrypto;



Ciphertext<DCRTPoly> replicate
(
    const Ciphertext<DCRTPoly> &input,
    const size_t size,
    const size_t initialPos,
    const size_t step=1
)
{
    Ciphertext<DCRTPoly> result = input;
    const size_t log2sizeF = std::floor(std::log2(size));
    size_t sizeLeft = size - (1 << log2sizeF);
    size_t adjustedPos = (initialPos < size / 2) ? initialPos
                                                 : initialPos - sizeLeft;
    Ciphertext<DCRTPoly> partial[log2sizeF + 1];
    partial[0] = result;
    for (size_t l = 0; l < log2sizeF; l++)
    {
        if ((adjustedPos >> l) & 1)
            result += result << step * (1 << l);
        else
            result += result >> step * (1 << l);
        partial[l + 1] = result;
    }
    while (sizeLeft > 0)
    {
        size_t highestPow2 = std::floor(std::log2(sizeLeft));
        if (initialPos < size / 2)
        {
            size_t adjustingIndex = step * (
                size - sizeLeft +
                (adjustedPos & ((1 << highestPow2) - 1)) - initialPos);
            result += partial[highestPow2] >> adjustingIndex;
        }
        else
        {
            size_t adjustingIndex = step * (
                initialPos - (adjustedPos & ((1 << highestPow2) - 1)) -
                sizeLeft + (1 << highestPow2));
            result += partial[highestPow2] << adjustingIndex;
        }
        sizeLeft -= 1 << highestPow2;
    }
    return result;
}


Ciphertext<DCRTPoly> sum
(
    const Ciphertext<DCRTPoly> &input,
    const size_t size,
    const size_t step=1
)
{
    Ciphertext<DCRTPoly> result = input;
    const size_t log2sizeF = std::floor(std::log2(size));
    Ciphertext<DCRTPoly> partial[log2sizeF + 1];
    partial[0] = result;
    for (size_t l = 0; l < log2sizeF; l++)
    {
        result += result << step * (1 << l);
        partial[l + 1] = result;
    }
    for (size_t l = 0; l < log2sizeF; l++)
        if ((size >> l) & 1)
            result += partial[l] << step * (size - (size & ((1 << (l + 1)) - 1)));
    return result;
}


int main(int argc, char* argv[])
{

    // Check if input and output file names are provided
    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0] << " <id> <verbose> <shared_memory_size> <serialized_size> <n_dim> <ind_cheby>\n";
        return 1;
    }

    const size_t id = std::stoul(argv[1]);
    const bool verbose = std::stoi(argv[2]);
    const size_t sharedMemorySize = std::stoul(argv[3]);
    const size_t serializedSize = std::stoul(argv[4]);
    const size_t ndimB = std::stoul(argv[5]);
    const bool indCheby = std::stoi(argv[6]);

    std::cout << "Child process " << id << ": START" << std::endl;


    ////////////////////////////////////////////////////////////////////////
    //                   Preparation (data-independent)                   //
    ////////////////////////////////////////////////////////////////////////
    
    // Load resources
    std::string storagePath = "data/tmp";
    const size_t n = Serial::deserializeFromFile<size_t>(storagePath + "/n");
    const size_t k = Serial::deserializeFromFile<size_t>(storagePath + "/k");
    CryptoContext<DCRTPoly> cryptoContext = Serial::deserializeFromFile<CryptoContext<DCRTPoly>>(storagePath + "/crypto-context");
    PublicKey<DCRTPoly> publicKey = Serial::deserializeFromFile<PublicKey<DCRTPoly>>(storagePath + "/public-key");
    // PrivateKey<DCRTPoly> secretKey = Serial::deserializeFromFile<PrivateKey<DCRTPoly>>(storagePath + "/secret-key");
    EvalKey<DCRTPoly> relinKey = Serial::deserializeFromFile<EvalKey<DCRTPoly>>(storagePath + "/relin-key");
    std::shared_ptr<std::map<usint, EvalKey<DCRTPoly>>> rotationKeys = Serial::deserializeFromFile<std::shared_ptr<std::map<usint, EvalKey<DCRTPoly>>>>(storagePath + "/rotation-keys");
    std::vector<double> cmpCoeffs = Serial::deserializeFromFile<std::vector<double>>(storagePath + "/cmp-coeffs");
    std::vector<double> indCoeffs = Serial::deserializeFromFile<std::vector<double>>(storagePath + "/ind-coeffs");
    
    cryptoContext->InsertEvalMultKey({relinKey});
    cryptoContext->InsertEvalAutomorphismKey(rotationKeys);
    const usint numSlots = cryptoContext->GetCryptoParameters()->GetEncodingParams()->GetBatchSize();
    const usint maxNumEl = MIN(numSlots / (k * k), n);
    const usint maxNumElCompressed = (numSlots / (k * k)) * (k * k);
    Ciphertext<DCRTPoly> zeroC = cryptoContext->Encrypt(
        publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(std::vector<double> {0.0})
    );
    size_t v = id;
    size_t u = v / (k * k); // index of the compressed ciphertext
    size_t j = v % (k * k);
    size_t j1 = j / k;
    size_t j2 = j % k;
    std::vector<double> spreadMask(maxNumElCompressed, 0.0);
    for (size_t w = 0; w < maxNumEl; w++)
        spreadMask[j1 * maxNumEl * k + w * k + j2] = 1.0;
    std::vector<double> blockMask(maxNumEl * k * k, 0.0);
    for (size_t j1 = 0; j1 < k; j1++)
        for (size_t j2 = 0; j2 < k; j2++)
            blockMask[j1 * maxNumEl * k + j2] = 1.0;
    Plaintext plaintext;
    std::vector<double> plaintextVec;

    // Map shared memory
    int shm_fd = shm_open("/shared_memory", O_RDWR, S_IRUSR | S_IWUSR);
    if (shm_fd == -1)
    {
        perror("shm_open");
        return 1;
    }
    char* sharedMem = static_cast<char*>(mmap(NULL, sharedMemorySize, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0));
    if (sharedMem == MAP_FAILED) {
        perror("mmap");
        return 1;
    }

    // Signal the main process that we are ready
    std::cout << "Child process " << id << ": READY FOR DATA" << std::endl;
    createFile(storagePath + "/proc-ready-for-data-" + std::to_string(id));



    ////////////////////////////////////////////////////////////////////////
    //                    Preparation (data-dependent)                    //
    ////////////////////////////////////////////////////////////////////////


    ////// WAIT FOR SIGNAL TO START

    while (!doesFileExist(storagePath + "/start-preparation"))
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    std::cout << "Child process " << id << ": PREPARING" << std::endl;


    ////// LOAD DATAPOINTS

    std::vector<std::vector<double>> xA = Serial::deserializeFromFile<std::vector<std::vector<double>>>(storagePath + "/xA");
    const size_t ndimA = xA[0].size();
    const size_t ndim = ndimA + ndimB;
    std::vector<Ciphertext<DCRTPoly>> XB(ndimB);
    for (size_t d = 0; d < ndimB; d++)
        XB[d] = Serial::deserializeFromFile<Ciphertext<DCRTPoly>>(storagePath + "/XB-" + std::to_string(u) + "-" + std::to_string(d));


    ////// ENCODE xA

    std::vector<std::vector<double>> xAreplR(ndimA);
    for (size_t d = 0; d < ndimA; d++)
    {
        xAreplR[d] = std::vector<double>(k * maxNumEl);
        for (size_t w = 0; w < maxNumEl; w++)
        {
            size_t i = v * maxNumEl + w;
            for (size_t j = 0; j < k; j++)
                xAreplR[d][w * k + j] = xA[i][d];
            if (i == n - 1) break;
        }
    }


    ////// RE-ENCODE XB

    Ciphertext<DCRTPoly> XBspread;
    Ciphertext<DCRTPoly> XBreplR[ndimB];
    Ciphertext<DCRTPoly> XBrepl[ndimB];
    for (size_t d = 0; d < ndimB; d++)
    {
        XBspread = XB[d] * spreadMask;
        XBreplR[d] = replicate(XBspread, k, j2);
        if (j1 > 0)
            XBreplR[d] = XBreplR[d] << (j1 * maxNumEl * k);
        XBrepl[d] = replicate(XBreplR[d], k, 0, maxNumEl * k);
    }

    // if (verbose)
    // {
    //     std::cout << std::fixed << std::setprecision(3);

    //     std::cout << "XBrepl : " << std::endl;
    //     std::cout << XBrepl[0]->GetLevel() << std::endl;
    //     cryptoContext->Decrypt(secretKey, XBrepl[0], &plaintext);
    //     plaintext->SetLength(numSlots);
    //     plaintextVec = plaintext->GetRealPackedValue();
    //     for (size_t j1 = 0; j1 < k; j1++)
    //     {
    //         for (size_t w = 0; w < maxNumEl; w++)
    //             for (size_t j2 = 0; j2 < k; j2++)
    //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;

    //     std::cout << "XBreplR : " << std::endl;
    //     std::cout << XBreplR[0]->GetLevel() << std::endl;
    //     cryptoContext->Decrypt(secretKey, XBreplR[0], &plaintext);
    //     plaintext->SetLength(numSlots);
    //     plaintextVec = plaintext->GetRealPackedValue();
    //     for (size_t j1 = 0; j1 < k; j1++)
    //     {
    //         for (size_t w = 0; w < maxNumEl; w++)
    //             for (size_t j2 = 0; j2 < k; j2++)
    //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }


    ////// SIGNAL THE MAIN PROCESS THAT WE ARE READY

    std::cout << "Child process " << id << ": READY" << std::endl;
    createFile(storagePath + "/proc-ready-" + std::to_string(id));



    ////////////////////////////////////////////////////////////////////////
    //                         LLoyd's Algorithm                          //
    ////////////////////////////////////////////////////////////////////////


    bool hasConverged = false;
    while (!hasConverged)
    {
        // Wait for trigger to start or to stop
        while (!doesFileExist(storagePath + "/start-computing"))
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if (doesFileExist(storagePath + "/stop"))
            {
                // Stop the child process
                std::cout << "Child process " << id << ": STOP" << std::endl;

                // Unmap shared memory and close the descriptor
                munmap(sharedMem, sharedMemorySize);
                close(shm_fd);

                return 0;
            }
        }

        // Get starting centers
        std::vector<std::vector<double>> cA = Serial::deserializeFromFile<std::vector<std::vector<double>>>(storagePath + "/cA");
        std::vector<std::vector<double>> cB = Serial::deserializeFromFile<std::vector<std::vector<double>>>(storagePath + "/cB");

        std::cout << "Child process " << id << ": COMPUTE" << std::endl;

        if (verbose)
        {
            std::cout << "Centers Alice: " << cA << std::endl;
            std::cout << "Centers Bob  : " << cB << std::endl << std::endl;
        }
        
        Ciphertext<DCRTPoly> sumA, sumB, clusterSize;

        // Alice computes her component of the distance of the datapoints from
        // each center and encodes it.
        std::vector<double> distAreplR(maxNumEl * k * k, 0.0);
        std::vector<double> distAreplC(maxNumEl * k * k, 0.0);
        for (size_t j1 = 0; j1 < k; j1++)
            for (size_t w = 0; w < maxNumEl; w++)
            {
                size_t i = v * maxNumEl + w;
                for (size_t j2 = 0; j2 < k; j2++)
                {
                    size_t index = j1 * maxNumEl * k + w * k + j2;
                    distAreplR[index] = 0.0;
                    distAreplC[index] = 0.0;
                    for (size_t d = 0; d < ndimA; d++)
                    {
                        distAreplR[index] += std::pow(cA[j2][d] - xA[i][d], 2);
                        distAreplC[index] += std::pow(cA[j1][d] - xA[i][d], 2);
                    }
                }
                if (i == n - 1) break;
            }

        // Alice encodes Bob's component of the centers and computes Bob's
        // component of the distance.
        Ciphertext<DCRTPoly> distBreplR = zeroC;
        Ciphertext<DCRTPoly> distBreplC = zeroC;
        for (size_t d = 0; d < ndimB; d++)
        {
            std::vector<double> cBreplR(maxNumEl * k * k, 0.0);
            std::vector<double> cBreplC(maxNumEl * k * k, 0.0);
            for (size_t j1 = 0; j1 < k; j1++)
                for (size_t w = 0; w < maxNumEl; w++)
                {
                    size_t i = v * maxNumEl + w;
                    for (size_t j2 = 0; j2 < k; j2++)
                    {
                        cBreplR[j1 * maxNumEl * k + w * k + j2] = -cB[j2][d];
                        cBreplC[j1 * maxNumEl * k + w * k + j2] = -cB[j1][d];
                    }
                    if (i == n - 1) break;
                }
            Ciphertext<DCRTPoly> tmpR = XBrepl[d] + cBreplR;
            Ciphertext<DCRTPoly> tmpC = XBrepl[d] + cBreplC;
            distBreplR += tmpR * tmpR;
            distBreplC += tmpC * tmpC;
        }

        // Alice aggregate the two components of the distance.
        Ciphertext<DCRTPoly> distR = (1.0 / ndim) * (distBreplR + distAreplR);
        Ciphertext<DCRTPoly> distC = (1.0 / ndim) * (distBreplC + distAreplC);

        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, distR, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "distR : " << std::endl;
        //     std::cout << distR->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }

        //     cryptoContext->Decrypt(secretKey, distC, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "distC : " << std::endl;
        //     std::cout << distC->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }
        
        // Alice finds the closest center to each datapoint (argmin of the
        // distances).
        auto start = std::chrono::high_resolution_clock::now();
        Ciphertext<DCRTPoly> cmpResult = cryptoContext->EvalChebyshevSeries(distR - distC, cmpCoeffs, -1.0, 1.0);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "cheby " << elapsed_seconds.count() << std::endl;

        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, cmpResult, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "cmpResult : " << std::endl;
        //     std::cout << cmpResult->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        // // Tie correction
        // std::vector<double> triangularMask(maxNumEl * k * k, 0.0);
        // for (size_t j1 = 0; j1 < k; j1++)
        //     for (size_t w = 0; w < maxNumEl; w++)
        //         for (size_t j2 = 0; j2 < k; j2++)
        //             if (j2 >= j1) triangularMask[j1 * maxNumEl * k + w * k + j2] = 1.0;
        // if (verbose)
        // {
        //     std::cout << "triangularMask : " << std::endl;
        //     std::cout << triangularMask->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << triangularMask[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        // Ciphertext<DCRTPoly> e = 4 * (1 - cmpResult) * cmpResult;
        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, e, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "e : " << std::endl;
        //     std::cout << e->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        // Ciphertext<DCRTPoly> correctionOffset1 = e;
        // for (size_t j = 1; j < k; j++)
        //     correctionOffset1 += e << j * maxNumEl * k;
        // e = e * triangularMask;
        // Ciphertext<DCRTPoly> correctionOffset2 = e;
        // for (size_t j = 1; j < k; j++)
        //     correctionOffset2 += e << j * maxNumEl * k;
        // Ciphertext<DCRTPoly> correctionOffset = correctionOffset2 - 0.5 * correctionOffset1;

        Ciphertext<DCRTPoly> rank = sum(cmpResult, k, maxNumEl * k);
        // rank += correctionOffset + (-0.5);

        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, rank, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     double max = *std::max_element(plaintextVec.begin(), plaintextVec.end());
        //     double min = *std::min_element(plaintextVec.begin(), plaintextVec.end());
        //     std::cout << "rank : " << std::endl;
        //     std::cout << "min: " << min << " - max: " << max << std::endl;
        //     std::cout << rank->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        Ciphertext<DCRTPoly> argmin = indCheby ?
            cryptoContext->EvalChebyshevSeries(rank, indCoeffs, -0.5, k) :
            cryptoContext->EvalPoly(rank, indCoeffs);

        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, argmin, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "argmin : " << std::endl;
        //     std::cout << argmin->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        std::vector<Ciphertext<DCRTPoly>> clusterA(ndimA);
        std::vector<Ciphertext<DCRTPoly>> clusterB(ndimB);
        for (size_t d = 0; d < ndimA; d++)
            clusterA[d] = argmin * xAreplR[d];
        for (size_t d = 0; d < ndimB; d++)
            clusterB[d] = argmin * XBreplR[d];
        clusterSize = argmin * std::vector<double>(k * maxNumEl, 1.0);

        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, clusterA[0], &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "clusterA : " << std::endl;
        //     std::cout << clusterA[0]->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }

        //     cryptoContext->Decrypt(secretKey, clusterB[0], &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "clusterB : " << std::endl;
        //     std::cout << clusterB[0]->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }

        //     cryptoContext->Decrypt(secretKey, clusterSize, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "clusterSize : " << std::endl;
        //     std::cout << clusterSize->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        Ciphertext<DCRTPoly> result = clusterSize
            + (clusterA[0] >> (maxNumEl * k))
            + (clusterB[0] >> (maxNumEl * 2 * k));
        result = sum(result, maxNumEl, k);
        for (size_t d = 1; d < MAX(ndimA, ndimB); d++)
        {
            if (d == 1)
                result = result * blockMask;
            Ciphertext<DCRTPoly> resultTmp = zeroC;
            if (d < ndimA)
                resultTmp += clusterA[d] >> (maxNumEl * k);
            if (d < ndimB)
                resultTmp += clusterB[d] >> (maxNumEl * 2 * k);
            resultTmp = sum(resultTmp, maxNumEl, k);
            result = result + ((resultTmp * blockMask) >> (d * k));
        }

        // if (verbose)
        // {
        //     cryptoContext->Decrypt(secretKey, result, &plaintext);
        //     plaintext->SetLength(numSlots);
        //     plaintextVec = plaintext->GetRealPackedValue();
        //     std::cout << "result : " << std::endl;
        //     std::cout << result->GetLevel() << std::endl;
        //     for (size_t j1 = 0; j1 < k; j1++)
        //     {
        //         for (size_t w = 0; w < maxNumEl; w++)
        //             for (size_t j2 = 0; j2 < k; j2++)
        //                 std::cout << plaintextVec[j1 * maxNumEl * k + w * k + j2] << " ";
        //         std::cout << std::endl;
        //     }
        // }

        // Write output
        const std::string result_s = Serial::serialize(result);
        assert(result_s.size() <= serializedSize);
        memcpy(
            sharedMem + serializedSize * id,
            result_s.c_str(),
            result_s.size()
        );

        // Signal that the output is ready
        createFile(storagePath + "/output-status-" + std::to_string(id));

        std::cout << "Child process " << id << ": DONE" << std::endl;
    }

    return 0;

}
