#ifndef INCLUDED_CRC32_H
#define INCLUDED_CRC32_H

#include <cstring>
#include <cstdio>
#include <stdint.h>

typedef uint32_t crc32_t;

template<typename T>
inline T reverse_crc_T(T data) {
  const int nBits = 8 * sizeof(T);
  T reflection = 0;
  for (int bit = 0; bit < nBits; ++bit) {
    if (data & 0x1) {
      reflection |= (1 << ((nBits - 1) - bit));
    }
    data = (data >> 1);
  }
  return reflection;
}

struct ReverseByteTable {
  unsigned char table[256];
  //
  ReverseByteTable() {
    init();
  }
  //
  void init() {
    for (int i = 0; i < 256; i++) {
      table[i] = reverse_crc_T((unsigned char)i);
    }
  }
};

inline unsigned char reverse_crc_byte(const unsigned char data) {
  const static ReverseByteTable reverseByteTable;
  return reverseByteTable.table[data];
}

struct CrcTable {
  crc32_t table[256];
  //
  CrcTable() {
    init();
  }
  //
  void init() {
    const crc32_t POLYNOMIAL = 0x04C11DB7;
    const int CRC_WIDTH = 8 * sizeof(crc32_t);
    const int CRC_TOPBIT = 1 << (CRC_WIDTH - 1);
    for (int dividend = 0; dividend < 256; ++dividend) {
      crc32_t remainder = dividend << (CRC_WIDTH - 8);
      for (unsigned char bit = 8; bit > 0; --bit) {
        if (remainder & CRC_TOPBIT) {
          remainder = (remainder << 1) ^ POLYNOMIAL;
        } else {
          remainder = (remainder << 1);
        }
      }
      table[dividend] = remainder;
    }
  }
};

inline crc32_t crc32(const crc32_t initial, const void* smessage, int nBytes) {
  // Compute CRC-32 of smessage
  // Code is originally from Michael Barr
  // http://www.barrgroup.com/Embedded-Systems/How-To/CRC-Calculation-C-Code
  const static CrcTable crcTable;
  const int CRC_WIDTH = 8 * sizeof(crc32_t);
  const unsigned char* message = (unsigned char*)smessage;
  crc32_t remainder = reverse_crc_T(initial) ^ 0xFFFFFFFF;
  for (int byte = 0; byte < nBytes; ++byte) {
    const unsigned char data = reverse_crc_byte(message[byte]) ^ (remainder >> (CRC_WIDTH - 8));
    remainder = crcTable.table[data] ^ (remainder << 8);
  }
  return reverse_crc_T(remainder) ^ 0xFFFFFFFF;
}

inline crc32_t crc32(const void* smessage, int nBytes) {
  return crc32(0, smessage, nBytes);
}

inline void crc32_check() {
  const char* test = "123456789";
  const crc32_t CHECK_VALUE = 0xCBF43926;
  std::printf("The check value for the %s standard is 0x%X\n", "CRC32", 0xCBF43926);
  std::printf("The crc32() of \"123456789\" is 0x%X\n", crc32(test, std::strlen(test)));
}

#endif
