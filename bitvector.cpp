#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cassert>
#include <array>
#include <iostream>
#include <boost/integer.hpp>

using namespace std;
using namespace boost;

typedef uint64_t u64;
#define U64_BIT (sizeof(u64) * CHAR_BIT)

#define LOGLOGN (4)

#define LOGN (1 << LOGLOGN)
#define LOG2N (LOGN * LOGN)
#define LOG4N (LOG2N * LOG2N)
#define HLOGN (LOGN / 2)
#define N (1 << LOGN)

typedef uint_t<LOGN>::exact r1_t;
typedef uint_t<2 * LOGLOGN>::exact r2_t;

template<unsigned size>
struct bitvector {
  array<u64, size / U64_BIT> data;
  unsigned operator[](u64 i) const {
    return (data[i / U64_BIT] & (1 << (i % U64_BIT))) >> (i % U64_BIT);
  }
  void set(u64 i) {
    data[i / U64_BIT] |= (1 << (i % U64_BIT));
  }
  unsigned popcount(u64 i, u64 j) const {
    assert(j - i < U64_BIT);
    assert(i / U64_BIT == j / U64_BIT);
    u64 t = data[i / U64_BIT];
    // cout << "i: " << i << " j: " << j << endl;
    // printf("t: %llx\n", t);
    i %= U64_BIT;
    j %= U64_BIT;
    t &= -(1 << i);
    // printf("t: %llx\n", t);
    t &= (1 << j) - 1;
    // printf("t: %llx\n", t);
    return __builtin_popcount(t);
  };
};

bitvector<N> B;
array<r1_t, N / LOG2N> R1;
array<r2_t, N / HLOGN> R2;

void init() {
  for (u64 i = 0; i < sizeof R1 / sizeof R1[0]; ++i) {
    if (i == 0) {
      R1[0] = 0;
      continue;
    }
    R1[i] = R1[i - 1];
    for (u64 j = 0; j < LOG2N; ++j) {
      R1[i] += B[(i - 1) * LOG2N + j];
    }
  }
  for (u64 i = 0; i < sizeof R2 / sizeof R2[0]; ++i) {
    if (i % (LOG2N / HLOGN) == 0) {
      R2[i] = 0;
      continue;
    }
    R2[i] = R2[i - 1];
    for (u64 j = 0; j < HLOGN; ++j) {
      R2[i] += B[(i - 1) * HLOGN + j];
    }
  }
}

u64 rank1(u64 x) {
  u64 r1 = R1[x/LOG2N];
  u64 r2 = R2[x/HLOGN];
  u64 r3 = B.popcount(x / HLOGN * HLOGN, x);
  // cout << r1 << " " << r2 << " " << r3 << endl;
  return r1 + r2 + r3;
}

u64 rank0(u64 x) {
  return x - rank1(x) + 1;      // 'x + 1 - rank(x)' can overflow
}

u64 select(u64 i) {
  u64 q = i / LOG2N;
  u64 r = i % LOG2N;
  if (S[q].m > LOG4N) {
    return S[q].u.T[r];
  } else {
    node *n = S[q].u.R;
    while (n->isbranch) {
      
    }
  }
}

int main()
{
  cout << "N: " << N << endl;

  B.set(10);
  B.set(15);
  B.set(20);
  B.set(25);
  B.set(10000);
  init();
  cout << rank1(7) << endl;
  cout << rank1(17) << endl;
  cout << rank1(27) << endl;
  cout << rank1(37) << endl;
  cout << rank1(10000) << endl;
  cout << rank1(10001) << endl;
}
