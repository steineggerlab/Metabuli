//
// Created by KJB on 26/08/2020.
//

#ifndef ADKMER3_MMAP_H
#define ADKMER3_MMAP_H
#include <cerrno>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>
#include <time.h>

template <typename T>
struct MmapedData
{
    T* data;
    size_t fileSize;
};

template <typename T>
MmapedData<T> mmapData(const char* filename, int mode = 1)
{
    struct MmapedData<T> mmapedData;
    struct stat stat1;
    int file;
    if(mode == 2) {
        file = open(filename, O_RDONLY);
    } else {
        file = open(filename, O_CREAT | O_RDWR);
    }
    int a;
    a = stat(filename, &stat1);
    mmapedData.fileSize = stat1.st_size;
    if(mode == 1) {
        mmapedData.data = static_cast<T *>(mmap(0, stat1.st_size + sizeof(T) * 2, PROT_WRITE | PROT_READ, MAP_SHARED,
                                                file, 0));
    } else if (mode == 2) { // Only read
        mmapedData.data = static_cast<T *>(mmap(0, stat1.st_size + sizeof(T) * 2, PROT_READ, MAP_SHARED,
                                                file, 0));
    }
    close(file);
    if(a == -1){
        mmapedData.fileSize = 0;
    }
    return mmapedData;
}


#endif //ADKMER3_MMAP_H
