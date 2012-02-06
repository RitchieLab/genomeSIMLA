
#include "snppool.h"
using namespace Genetics;
int main() {
	SnpPool::Initialize(5, 5, 4);
	SnpPool *pool = SnpPool::Instance();
	SnpAligned *snp = pool->GetSnp(4);
	pool->ReleaseSnp(snp);
	pool->Release();
	return 0;
}
