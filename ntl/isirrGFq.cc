//  g++ -O3 -o isirrGFq isirrGFq.cc -l ntl

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pEX.h>

using namespace std;
using namespace NTL;

int main()
{
   ZZ p;
   cin >> p;
   ZZ_p::init(p);

   ZZ_pX f;
   cin >> f;
   ZZ_pE::init(f);

   ZZ_pEX h;
   cin >> h;

   if (IterIrredTest(h)) 
     cout << "true\n";
   else
     cout << "false\n";
     
}
