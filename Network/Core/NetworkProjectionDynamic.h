#include "Network.h"

/// <summary>	Not used. </summary>

class ProjectionDynamic : public Projection
{
public:

	ProjectionDynamic();

	void SimulateEvent(UnitModifier* e);
	void ModifyProjection();

	void AddProjection(Unit* pre, Unit* post, bool firstRun);
	
protected:

	double tau_psc_;   //!< [ms] time constant of postsyn current
	double tau_fac_;   //!< [ms] time constant for fascilitation
	double tau_rec_;   //!< [ms] time constant for recovery
	double U_;         //!< asymptotic value of probability of release
	double x_;         //!< amount of resources in recovered state
	double y_;         //!< amount of resources in active state
	double u_;         //!< actual probability of release	
};