#include "PhyElementFrame.h"
#include "PhyMaterial.h"
#include "PhyNode.h"

void PhyElementFrame::setInternalMaterialProperties(PhyMaterial* pMat)
{
	E = pMat->paras(mpb_E);
	A = pMat->paras(mpb_A);
	I = pMat->paras(mpb_I);
//	eType = 4 Frame matfeng = [E A I] (454d)
}

void PhyElementFrame::setGeometry()
{
	// Complete
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


void PhyElementFrame::Calculate_ElementStiffness_Force()
{
	a1 = E * A / L;
	a2 = E * I / pow(L, 3);

	//! 1. stiffness matrix in local coordinate system
	// kbar
	kLocalCoordinate.resize(6, 6);
	kLocalCoordinate = 0.0;
	kLocalCoordinate(0, 0) = a1;
	kLocalCoordinate(0, 3) = -a1;
	kLocalCoordinate(1, 1) = 12.0 * a2;
	kLocalCoordinate(1, 2) = 6 * L * a2;
	kLocalCoordinate(1, 4) = -12.0 * a2;
	kLocalCoordinate(1, 5) = 6 * L * a2;
	kLocalCoordinate(2, 1) = 6 * L * a2;
	kLocalCoordinate(2, 2) = 4 * pow(L,2) * a2;
	kLocalCoordinate(2, 4) = -6 * L * a2;
	kLocalCoordinate(2, 5) = 2 * pow(L,2) * a2;
	kLocalCoordinate(3, 0) = -a1;
	kLocalCoordinate(3, 3) = a1;
	kLocalCoordinate(4, 1) = -12.0 * a2;
	kLocalCoordinate(4, 2) = -6 * L * a2;
	kLocalCoordinate(4, 3) = 12.0 * a2;
	kLocalCoordinate(4, 4) = -6 * L * a2;
	kLocalCoordinate(5, 1) = 6 * L * a2;
	kLocalCoordinate(5, 2) = 2 * pow(L,2) * a2;
	kLocalCoordinate(5, 4) = -6 * L * a2;
	kLocalCoordinate(5, 5) = 4 * pow(L,2) * a2;

	//! 2. Transformation matrix T
	T.resize(6, 6);
	T(0, 0) = c;
	T(0, 1) = s;
	T(1, 0) = -s;
	T(1, 1) = c;
	T(0, 0) = c;
	T(0, 1) = s;
	T(2, 2) = 1;

	T(3, 3) = c;
	T(3, 4) = s;
	T(4, 3) = -s;
	T(4, 4) = c;
	T(3, 3) = c;
	T(3, 4) = s;
	T(5, 5) = 1;

	//step 3: calculate k as k=T transpose * kbar  * T
	/*
	T_transpose.resize(6, 6);
	for(int xi=0;xi<6;i++){
		for(int xj=0;xj<6;j++)
		{
				T_transpose[j][i]=T[i][j];
		}
	 }
	 */

	//ke = T_transpose * kLocalCoordinate* T;
	/*
	MATRIX	mult.resize(6,6)
	MATRIX	ke.resize(6,6)
		for(bi = 0; bi < T_transpose.rows(); ++bi)
			for(bj = 0; bj < kLocalCoordinate.columns(); ++bj)
					for(bk = 0; bk < T_transpose.columns(); ++bk)
					{
							mult[bi][bj] += a[bi][bk] * b[bk][bj];
					}

	for(ai = 0; ai < mult.rows(); ++ai)
			for(aj = 0; aj < T.columns(); ++aj)
					for(ak = 0; ak < mult.columns(); ++ak)
					{
							ke[bi][bj] += a[bi][bk] * b[bk][bj];
					}
					*/
					ke.resize(6,6);
					ke = 0.0;
					for (int i = 0; i < 6; ++i)
						for (int j = 0; i < 6; ++i)
							for (int k = 0; i < 6; ++i)
								for (int l = 0; i < 6; ++i)
									ke(i,j) += T(k,i) * kLocalCoordinate(k,1) * T(l,j);


}


void PhyElementFrame::SpecificOutput(ostream& out) const
{
	// Complete
}
