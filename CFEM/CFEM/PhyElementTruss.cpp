#include "PhyElementTruss.h"
#include "PhyMaterial.h"
#include "PhyNode.h"

void PhyElementTruss::setInternalMaterialProperties(PhyMaterial* pMat)
{
A = pMat -> paras(mpb_A);
E = pMat -> paras(mpb_E);
}

void PhyElementTruss::setGeometry()
{
	VECTOR *crd0, *crd1;
	crd1 = &eNodePtrs[1]->coordinate;
	crd0 = &eNodePtrs[0]->coordinate;

	int sz = crd1->size();
	if (sz != 2)
		THROW("implementation only for 2D truss");
	double delX, delY;
	delX = (*crd1)(0) - (*crd0)(0);
	delY = (*crd1)(1) - (*crd0)(1);
	L = sqrt(delX * delX + delY * delY);
	c = delX / L;
	s = delY / L;
}


void PhyElementTruss::Calculate_ElementStiffness_Force()
{
	//step 3: calculate k as k=T transpose * kbar  * T
	ke.resize(4,4);

	double factor = A * E /L;
	//for (int I = 0; I<2; ++I)
		//for (int J = 0; J<2; ++J)
		//{
			//double f2 = factor;
			//if ((I + J) % 2 != 0)
				//f2 = -factor;
			//ke(I,J) = c*c*f2;
			//ke(I+1,J) = ke(I,J+1) = c*s*f2;
			//ke(I+1,J+1) = s*s*f2;
			ke(0,0) =factor*c*c;
			ke(0,1) = factor*c*s;
			ke(0,2) =-factor*c*c;
			ke(0,3) =-factor*c*s;

			ke(1,0) =factor*c*s;
			ke(1,1) =factor*s*s;
			ke(1,2) =-factor*c*s;
			ke(1,3) =-factor*s*s;

			ke(2,0) =-factor*c*c;
			ke(2,1) =-factor*c*s;
			ke(2,2) =factor*c*c;
			ke(2,3) =factor*c*s;

			ke(3,0) =-factor*c*s;
			ke(3,1) =-factor*s*s;
			ke(3,2) = factor*c*s;
			ke(3,3) = factor*s*s;
			cout << "stiffness matrix is" <<  ke << "\n";
		//}


}

void PhyElementTruss::SpecificOutput(ostream& out) const
{
	double T = A * E / L * (c * (edofs(2) - edofs(0)) + s * (edofs(3) - edofs(1)));
	out << T;
}
