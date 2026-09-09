#ifndef METABULI_INFOINDEX_H
#define METABULI_INFOINDEX_H

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "common.h"

// Metadata stored in db.parameters for the final "info" file. Older databases
// do not have these keys; absence means the file is the original uint32 stream.
struct InfoIndexMetadata {
    std::string format = "uint32"; // "uint32" or "packed_uint"
    uint8_t idBits = 32;           // Bits used by each logical ID.
    size_t idCount = 0;            // Logical ID count, excluding padded lanes.

    bool isPacked() const {
        return format == "packed_uint" && idBits < 32;
    }
};

struct InfoIndex {
    static constexpr uint8_t WORD_BITS = 64;

    // Metadata lives next to the final info file in db.parameters.
    static std::string dbParameterFile(const std::string &infoFileName) {
        const size_t slash = infoFileName.find_last_of('/');
        if (slash == std::string::npos) {
            return "db.parameters";
        }
        return infoFileName.substr(0, slash + 1) + "db.parameters";
    }

    // Temporary merge chunks are named "<n>_info" and remain raw uint32. Only
    // the final database file is named exactly "info" and can use metadata.
    static bool isFinalInfoFile(const std::string &infoFileName) {
        const size_t slash = infoFileName.find_last_of('/');
        const std::string base = (slash == std::string::npos)
                                     ? infoFileName
                                     : infoFileName.substr(slash + 1);
        return base == "info";
    }

    static uint8_t idsPerWord(uint8_t idBits) {
        return std::max<uint8_t>(1, WORD_BITS / idBits);
    }

    static uint64_t maskForBits(uint8_t idBits) {
        return idBits == 64 ? UINT64_MAX : ((uint64_t{1} << idBits) - 1);
    }

    static size_t packedWordCount(size_t idCount, uint8_t idBits) {
        const size_t perWord = idsPerWord(idBits);
        return (idCount + perWord - 1) / perWord;
    }

    static uint8_t chooseIdBits(uint32_t maxId) {
        // Tiers are chosen only where they improve IDs-per-word density. For
        // example, 18 and 21 bits both fit three IDs per uint64, so 21 keeps
        // more future headroom at the same file size.
        static constexpr uint8_t TIERS[] = {8, 10, 12, 16, 21, 32};
        for (uint8_t bits : TIERS) {
            if (bits == 32 || maxId <= maskForBits(bits)) {
                return bits;
            }
        }
        return 32;
    }

    static InfoIndexMetadata makeMetadata(uint32_t maxId, size_t idCount) {
        InfoIndexMetadata meta;
        meta.idBits = chooseIdBits(maxId);
        meta.idCount = idCount;
        meta.format = (meta.idBits < 32) ? "packed_uint" : "uint32";
        return meta;
    }

    static InfoIndexMetadata loadMetadata(const std::string &infoFileName) {
        InfoIndexMetadata meta;
        // Non-final info files are intermediate merge inputs. They must stay
        // readable without a db.parameters side channel.
        if (!isFinalInfoFile(infoFileName)) {
            return meta;
        }

        std::ifstream in(dbParameterFile(infoFileName));
        if (!in.is_open()) {
            return meta;
        }

        std::string line;
        while (std::getline(in, line)) {
            const size_t tab = line.find('\t');
            if (tab == std::string::npos) {
                continue;
            }
            const std::string key = line.substr(0, tab);
            const std::string value = line.substr(tab + 1);
            if (key == "Info_format") {
                meta.format = value;
            } else if (key == "Info_id_bits") {
                meta.idBits = static_cast<uint8_t>(std::stoul(value));
            } else if (key == "Info_id_count") {
                meta.idCount = static_cast<size_t>(std::stoull(value));
            }
        }

        // Be conservative with incomplete or unknown metadata. Falling back to
        // uint32 preserves compatibility with old databases and temp chunks.
        if (meta.format != "packed_uint" || meta.idBits >= 32 || meta.idCount == 0) {
            meta.format = "uint32";
            meta.idBits = 32;
        }
        return meta;
    }

    static void appendMetadata(const std::string &parameterFileName,
                               const InfoIndexMetadata &meta) {
        FILE *handle = fopen(parameterFileName.c_str(), "a");
        if (handle == nullptr) {
            std::cerr << "Could not open " << parameterFileName
                      << " for writing info index metadata\n";
            std::exit(EXIT_FAILURE);
        }
        writeMetadata(handle, meta);
        fclose(handle);
    }

    static void upsertMetadata(const std::string &parameterFileName,
                               const InfoIndexMetadata &meta) {
        std::vector<std::string> keptLines;
        {
            std::ifstream in(parameterFileName);
            std::string line;
            while (std::getline(in, line)) {
                const size_t tab = line.find('\t');
                const std::string key = (tab == std::string::npos) ? line : line.substr(0, tab);
                // Maintenance commands may run more than once. Drop all old info
                // metadata keys and write one authoritative block at the end.
                if (!isMetadataKey(key)) {
                    keptLines.push_back(line);
                }
            }
        }

        const std::string tmpFileName = parameterFileName + ".info.tmp";
        FILE *handle = fopen(tmpFileName.c_str(), "w");
        if (handle == nullptr) {
            std::cerr << "Could not open " << tmpFileName
                      << " for writing info index metadata\n";
            std::exit(EXIT_FAILURE);
        }
        for (const std::string &line : keptLines) {
            fprintf(handle, "%s\n", line.c_str());
        }
        writeMetadata(handle, meta);
        fclose(handle);

        if (std::rename(tmpFileName.c_str(), parameterFileName.c_str()) != 0) {
            std::cerr << "Could not replace " << parameterFileName
                      << " with updated info index metadata\n";
            std::exit(EXIT_FAILURE);
        }
    }

private:
    static bool isMetadataKey(const std::string &key) {
        return key == "Info_format" ||
               key == "Info_id_bits" ||
               key == "Info_id_count";
    }

    static void writeMetadata(FILE *handle, const InfoIndexMetadata &meta) {
        fprintf(handle, "Info_format\t%s\n", meta.format.c_str());
        fprintf(handle, "Info_id_bits\t%u\n", static_cast<unsigned int>(meta.idBits));
        fprintf(handle, "Info_id_count\t%lu\n", static_cast<unsigned long>(meta.idCount));
    }
};

class InfoIndexWriter {
private:
    const uint8_t idBits;
    const uint8_t valuesPerWord;
    const uint64_t idMask;
    WriteBuffer<uint64_t> wordBuffer;
    uint64_t currentWord = 0; // Accumulates IDs until one uint64 word is full.
    uint8_t lane = 0;         // Logical slot inside currentWord.
    bool closed = false;

public:
    InfoIndexWriter(const std::string &fileName,
                    uint8_t idBits,
                    size_t wordBufferSize)
        : idBits(idBits),
          valuesPerWord(InfoIndex::idsPerWord(idBits)),
          idMask(InfoIndex::maskForBits(idBits)),
          wordBuffer(fileName, wordBufferSize) {}

    ~InfoIndexWriter() {
        close();
    }

    inline void write(uint32_t id) {
        // Validate once on write so the reader can stay branch-light.
        if (unlikely((static_cast<uint64_t>(id) & ~idMask) != 0)) {
            std::cerr << "Info ID " << id << " exceeds " << static_cast<int>(idBits)
                      << " packed bits\n";
            std::exit(EXIT_FAILURE);
        }

        // Hot packing path: one masked shift/or per logical ID, then flush a
        // complete 64-bit word when the current lane group is full.
        currentWord |= (static_cast<uint64_t>(id) << (lane * idBits));
        ++lane;
        if (unlikely(lane == valuesPerWord)) {
            wordBuffer.write(&currentWord);
            currentWord = 0;
            lane = 0;
        }
    }

    void close() {
        if (closed) {
            return;
        }
        // The last word may contain unused lanes; Info_id_count tells readers
        // where the logical stream ends.
        if (lane != 0) {
            wordBuffer.write(&currentWord);
            currentWord = 0;
            lane = 0;
        }
        wordBuffer.close();
        closed = true;
    }
};

class InfoIndexReader {
private:
    InfoIndexMetadata metadata;
    // Exactly one buffer is active. Keeping separate typed buffers avoids
    // per-ID virtual dispatch or conversions in DeltaIdxReader's hot loop.
    std::unique_ptr<ReadBuffer<uint32_t>> plainBuffer;
    std::unique_ptr<ReadBuffer<uint64_t>> packedBuffer;
    uint8_t valuesPerWord = 2;
    uint64_t idMask = UINT32_MAX;
    uint64_t currentWord = 0; // Current packed word being decoded.
    uint8_t lane = 0;         // Next logical ID slot inside currentWord.
    size_t readCount = 0;     // Logical IDs returned so far from this reader.
    size_t totalValueNum = 0;

    inline bool loadPackedWord() {
        // ReadBuffer batches disk IO; this only advances to the next uint64
        // word when all packed lanes from the current word have been consumed.
        if (unlikely(packedBuffer->p >= packedBuffer->end)) {
            if (packedBuffer->loadBuffer() == 0) {
                return false;
            }
        }
        currentWord = *packedBuffer->p++;
        lane = 0;
        return true;
    }

public:
    InfoIndexReader(const std::string &infoFileName, size_t readBufferSize)
        : metadata(InfoIndex::loadMetadata(infoFileName)) {
        if (metadata.isPacked()) {
            valuesPerWord = InfoIndex::idsPerWord(metadata.idBits);
            idMask = InfoIndex::maskForBits(metadata.idBits);
            // readBufferSize is historically a uint32 count. Halving keeps the
            // packed uint64 byte budget close to the old reader budget.
            const size_t wordBufferSize = std::max<size_t>(1, readBufferSize / 2);
            packedBuffer = std::make_unique<ReadBuffer<uint64_t>>(infoFileName,
                                                                  wordBufferSize);
            lane = valuesPerWord;
            totalValueNum = metadata.idCount;
        } else {
            plainBuffer = std::make_unique<ReadBuffer<uint32_t>>(infoFileName,
                                                                 readBufferSize);
            totalValueNum = FileUtil::getFileSize(infoFileName) / sizeof(uint32_t);
        }
    }

    bool isPacked() const {
        return metadata.isPacked();
    }

    size_t getTotalValueNum() const {
        return totalValueNum;
    }

    bool loadAt(size_t logicalOffset) {
        readCount = logicalOffset;
        if (!metadata.isPacked()) {
            return plainBuffer->loadBufferAt(logicalOffset) > 0 || logicalOffset == totalValueNum;
        }

        // Split files store logical ID offsets, not physical byte/word offsets.
        // Convert once during seek; next() then only shifts and masks.
        const size_t wordOffset = logicalOffset / valuesPerWord;
        const uint8_t targetLane = logicalOffset % valuesPerWord;
        if (packedBuffer->loadBufferAt(wordOffset) == 0 && logicalOffset != totalValueNum) {
            return false;
        }

        lane = valuesPerWord;
        currentWord = 0;
        if (targetLane != 0) {
            if (!loadPackedWord()) {
                return false;
            }
            lane = targetLane;
        }
        return true;
    }

    inline bool next(uint32_t &id) {
        if (!metadata.isPacked()) {
            if (unlikely(plainBuffer->p >= plainBuffer->end)) {
                if (plainBuffer->loadBuffer() == 0) {
                    return false;
                }
            }
            id = *plainBuffer->p++;
            ++readCount;
            return true;
        }

        if (unlikely(readCount >= totalValueNum)) {
            return false;
        }
        if (unlikely(lane >= valuesPerWord && !loadPackedWord())) {
            return false;
        }

        // Hot unpacking path: the split seek has already selected the word and
        // lane, so retrieval is one shift, one mask, and one lane increment.
        id = static_cast<uint32_t>((currentWord >> (lane * metadata.idBits)) & idMask);
        ++lane;
        ++readCount;
        return true;
    }
};

inline void packInfoFileInPlace(const std::string &infoFileName,
                                const InfoIndexMetadata &metadata) {
    if (!metadata.isPacked()) {
        return;
    }

    // Final writers first produce the old uint32 stream. Packing afterwards
    // lets merge code compute max ID/count in one pass without buffering all IDs.
    const std::string packedFileName = infoFileName + ".packed.tmp";
    {
        ReadBuffer<uint32_t> reader(infoFileName, 1024 * 1024 * 16);
        InfoIndexWriter writer(packedFileName, metadata.idBits, 1024 * 1024 * 8);
        uint32_t id = 0;
        size_t written = 0;
        while (reader.p < reader.end || reader.loadBuffer() != 0) {
            while (reader.p < reader.end) {
                id = *reader.p++;
                writer.write(id);
                ++written;
            }
        }
        // This catches accidental count mismatches before replacing the final
        // database info file.
        if (written != metadata.idCount) {
            std::cerr << "Info ID count changed while packing " << infoFileName
                      << ": expected " << metadata.idCount << ", saw " << written << "\n";
            std::exit(EXIT_FAILURE);
        }
        writer.close();
    }

    if (std::remove(infoFileName.c_str()) != 0 ||
        std::rename(packedFileName.c_str(), infoFileName.c_str()) != 0) {
        std::cerr << "Could not replace " << infoFileName
                  << " with packed info index\n";
        std::exit(EXIT_FAILURE);
    }
}

#endif // METABULI_INFOINDEX_H
