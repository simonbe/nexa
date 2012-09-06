#include <math.h>
#include "NetworkProjectionModifier.h"
#include "NetworkUnits.h"
//#include "NetworkUnitModifier.h"



void ProjectionModifier::Initialize(Projection* Projection)
{
	m_projection = Projection;

	for(int i = 0;i<m_parentPopulationModifier.size();i++)
	{
		m_parentPopulationModifier[i]->AddChildProjectionModifier(this);
	}

	for(int i=0;i<m_parentProjectionModifier.size();i++)
	{
		m_parentProjectionModifier[i]->AddChildProjectionModifier(this);
	}
}

UnitModifier* ProjectionModifier::GetTransferFunction()
{
	return m_transferFunction;
}

ProjectionModifier::~ProjectionModifier()
{
	if(m_transferFunction!=NULL)
		delete m_transferFunction;
	
/*	for(int i=0;i<m_parentProjectionModifier.size();i++)
		delete m_parentProjectionModifier[i];

	for(int i=0;i<m_childProjectionModifier.size();i++)
		delete m_childProjectionModifier[i];

	for(int i=0;i<m_parentPopulationModifier.size();i++)
		delete m_parentPopulationModifier[i];
		*/
}