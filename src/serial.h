// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

#include <sys/stat.h> // Required for struct stat


inline bool doesFileExist(
    const std::string& name
)
{
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}


void createFile(
    const std::string& filepath
)
{
    std::ofstream output(filepath);
}


void deleteFile(
    const std::string& filepath
)
{
    std::remove(filepath.c_str());
}


namespace lbcrypto::Serial {

    template <typename T>
    std::string serialize(
        const T& obj
    )
    {
        std::stringstream objSerialized;
        Serial::Serialize(obj, objSerialized, SerType::BINARY);
        return objSerialized.str();
    }

    template <typename T>
    T deserialize(
        const std::string objSerialized
    )
    {
        T obj;
        std::stringstream ss(objSerialized);
        Serial::Deserialize(obj, ss, SerType::BINARY);
        return obj;
    }

    template <typename T>
    void serializeToFile(
        const T& obj, 
        const std::string& filename
    )
    {
        std::ofstream outFile(filename, std::ios::binary);
        if (!outFile)
            throw std::runtime_error("Error opening file for writing: " + filename);
        Serial::Serialize(obj, outFile, SerType::BINARY);
        outFile.close();
    }

    template <typename T>
    T deserializeFromFile(
        const std::string& filename
    )
    {
        T obj;
        std::ifstream inFile(filename, std::ios::binary);
        if (!inFile)
            throw std::runtime_error("Error opening file for reading: " + filename);
        Serial::Deserialize(obj, inFile, SerType::BINARY);
        inFile.close();
        return obj;
    }

}