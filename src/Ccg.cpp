#include "BoussinesqFft.h"
#include "Ccg.h"

void ccg(sub, topo, offset, gradient, gap) {
  auto nodDeflections = BoussinesqFft(1.0, sub, topo);
  double delta = 0;
  double G_old = 1.0;
  auto gap = nodDeflections - topo;
  auto contactArea = getNotNullNodes();
  auto g_mean = contactArea * sumNotNullNodes();
}
