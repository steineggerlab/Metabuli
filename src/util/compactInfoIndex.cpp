#include "Command.h"
#include "FileUtil.h"
#include "InfoIndex.h"
#include "LocalParameters.h"
#include "validateDatabase.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

struct PlainInfoStats {
    uint32_t maxId = 0;
    size_t idCount = 0;
};

PlainInfoStats scanPlainInfoFile(const std::string &infoFileName) {
    const uint64_t fileSize = FileUtil::getFileSize(infoFileName);
    if (fileSize == 0) {
        std::cerr << "Error: info file is empty: " << infoFileName << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (fileSize % sizeof(uint32_t) != 0) {
        std::cerr << "Error: raw info file size is not a multiple of "
                  << sizeof(uint32_t) << ": " << infoFileName << std::endl;
        std::exit(EXIT_FAILURE);
    }

    PlainInfoStats stats;
    ReadBuffer<uint32_t> reader(infoFileName, 1024 * 1024 * 16);
    while (reader.p < reader.end || reader.loadBuffer() != 0) {
        while (reader.p < reader.end) {
            const uint32_t id = *reader.p++;
            stats.maxId = std::max(stats.maxId, id);
            ++stats.idCount;
        }
    }
    return stats;
}

int compactInfoIndexDatabase(const std::string &dbDir) {
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        std::cerr << "Error: database directory does not exist: " << dbDir << std::endl;
        return 1;
    }

    const std::string infoFileName = dbDir + "/info";
    const std::string parameterFileName = dbDir + "/db.parameters";
    if (!FileUtil::fileExists(infoFileName.c_str())) {
        std::cerr << "Error: info file is missing in database directory: "
                  << dbDir << std::endl;
        return 1;
    }

    const uint64_t originalBytes = FileUtil::getFileSize(infoFileName);
    const InfoIndexMetadata existingMetadata = InfoIndex::loadMetadata(infoFileName);
    if (existingMetadata.isPacked()) {
        const size_t expectedBytes =
            InfoIndex::packedWordCount(existingMetadata.idCount, existingMetadata.idBits) *
            sizeof(uint64_t);
        // If metadata and file size disagree, guessing the current physical
        // layout could corrupt the DB. Stop and let the user inspect it.
        if (originalBytes != expectedBytes) {
            std::cerr << "Error: db.parameters says info is packed, but info size is "
                      << originalBytes << " bytes; expected " << expectedBytes
                      << " bytes." << std::endl;
            return 1;
        }
        std::cout << "Info index is already packed: "
                  << static_cast<int>(existingMetadata.idBits) << " bits, "
                  << existingMetadata.idCount << " IDs." << std::endl;
        return validateDatabase(dbDir);
    }

    std::cout << "Scanning raw info index..." << std::endl;
    const PlainInfoStats stats = scanPlainInfoFile(infoFileName);
    InfoIndexMetadata metadata = InfoIndex::makeMetadata(stats.maxId, stats.idCount);
    std::cout << "Info ID count       : " << metadata.idCount << std::endl;
    std::cout << "Max info ID         : " << stats.maxId << std::endl;
    std::cout << "Selected info format: " << metadata.format << " ("
              << static_cast<int>(metadata.idBits) << " bits)" << std::endl;

    if (!metadata.isPacked()) {
        // IDs that need 32 bits cannot be compacted by this format. Still write
        // metadata so future readers know the logical count explicitly.
        InfoIndex::upsertMetadata(parameterFileName, metadata);
        std::cout << "Info IDs require 32 bits; info file was left unchanged." << std::endl;
        return validateDatabase(dbDir);
    }

    const size_t expectedPackedBytes =
        InfoIndex::packedWordCount(metadata.idCount, metadata.idBits) * sizeof(uint64_t);
    std::cout << "Original info bytes : " << originalBytes << std::endl;
    std::cout << "Packed info bytes   : " << expectedPackedBytes << std::endl;

    // packInfoFileInPlace streams the raw uint32 file, validates the logical ID
    // count, and only then replaces DBDIR/info with the packed file.
    packInfoFileInPlace(infoFileName, metadata);
    InfoIndex::upsertMetadata(parameterFileName, metadata);

    const uint64_t packedBytes = FileUtil::getFileSize(infoFileName);
    if (packedBytes != expectedPackedBytes) {
        std::cerr << "Error: packed info size is " << packedBytes
                  << " bytes; expected " << expectedPackedBytes << " bytes." << std::endl;
        return 1;
    }

    std::cout << "Info index compacted successfully." << std::endl;
    std::cout << "Saved bytes         : " << (originalBytes - packedBytes) << std::endl;
    return validateDatabase(dbDir);
}

} // namespace

int compactInfoIndex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    return compactInfoIndexDatabase(par.filenames[0]);
}
