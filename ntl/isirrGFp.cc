//  g++ -O3 -o isirrGFp isirrGFp.cc -l ntl
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

using namespace std;
using namespace NTL;

int main()
{
   ZZ p;
   cin >> p;
   ZZ_p::init(p);

   ZZ_pX f;
   GF2X g;

   if (p == 2) {
     cin >> g;
     if (IterIrredTest(g)) 
       cout << "true\n";
     else
       cout << "false\n";
   } else {
     cin >> f;
     if (IterIrredTest(f)) 
       cout << "true\n";
     else
       cout << "false\n";
   }
}
