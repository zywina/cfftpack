
#include <ql/math/primenumbers.hpp>
using namespace QuantLib;

#include <stdio.h>
#include <stdint.h>
#include <iostream>
using namespace std;



int main(){
  int N = 512;
  printf("const int NPRIMES = %d\n",N);
  printf("const unsinged int Primes[NPRIMES] = {\n\t");
  for (int x=0; x<N; x++){
    int y = PrimeNumbers::get(x);
    printf("%d, ",y);
    if (x>0 && x%20==0) printf("\n\t");
    if (x > 0xffffffffULL) break;
  }
  printf("};\n\n");
  return 0;
}
