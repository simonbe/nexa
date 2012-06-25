#include "NetworkUnitsIF.h"
#include "NetworkUnits.h"
#include <iostream>
#include <cmath>
#include <deque>


UnitIF::UnitIF()
{
	//m_unitType = "minicolumn";
	Km = 0.0041;
	Vth = 60000;
	Vahp = 0;
	Vk = 0;
	Ik = 0;
}

void UnitIF::SimulateEventQueue()
{
//	if(m_eventsIncoming.size()>0)
//	{
		//vector<long> currentIncomingSpikes = m_eventsIncoming[0];
		vector<long> currentIncomingSpikes;// = m_eventsIncoming[0];

		//if(currentIncomingSpikes.size()>0)
		//{
			float totWeights = 0;

			//			m_eventsIncoming[0].clear();

			vector<Connection*> conns = this->GetPopulation()->GetIncomingConnections();
			vector<UnitModifier*> eus;
			vector<float> weights;

			// put connection vars in network instead of connection to avoid check

			for(int m=0;m<conns.size();m++)
			{
				//for(int j=0;j<1;j++)//conns.size();j++)
				//{
				if(conns[m]->IsOn())
				{

					// currently fix after eventsIncoming removed, used as in graded units, even though incoming non-spikes should not need to be checked
					// change into storing incoming spike-ids in bufferSpikes or something similar

					vector<long> preIds = conns[m]->GetPreIds(this->GetUnitId());

					if(preIds.size()>0)
					{
						for(int j=0;j<preIds.size();j++)
						{

							//for(int i=0;i<currentIncomingSpikes.size();i++)//m_eventsIncoming.size();i++)
							//{
							float incomingBufferData = m_network->GetPreValue(preIds[j]);//GetIncomingBufferData(preIds[j]);//it->first);

							if(incomingBufferData!=0)
							{
								//UnitModifier* eu = network()->GetUnitModifierIncoming(currentIncomingSpikes[i]);

								//m_value += m_eventsIncoming[i]->GetValue() * conns[j]->GetWeight(m_eventsIncoming[i]->GetFromUnitId(),this->GetUnitId());
								bool exists = true;

								if(conns.size()>1) // if more than one incoming connection we need to know which connection it is coming from
								{
									exists = true;//network()->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());//conns[j]->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());
								}

								if(exists)
								{
									//totWeights += conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
									float weight = network()->GetWeight(preIds[j],this->GetUnitId());//eu->GetFromUnitId(),this->GetUnitId());//conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
									//if(weight<0) // inhibitory
									//{
									Ik+=weight*incomingBufferData;//weight*eu->GetEventData()[0];
									//if(eu->GetEventTypeId()==2)
//										bool b=false;
									//	S_.y_[DG_IN] += -weight * V_.g0_in_;
									//}
									//else // excitatory
									//	;//S_.y_[DG_EX] += weight * V_.g0_ex_;

									//weights.push_back(weight);
									//eus.push_back(eu);
									//eu->SetValue(eu->GetValue() * weight);
								}
							}
						}
					}
				}
			}

			//if(S_.y_[DG_EX] > 300 * V_.g0_ex_)
			//	S_.y_[DG_EX] = 300 * V_.g0_ex_;
			//}

			// elsewhere now
			//m_eventsIncoming.erase(m_eventsIncoming.begin(),m_eventsIncoming.begin()+1);
			//	}

	// matlab code
	/*	% output: spikes
	% Vk1: updated potential
	% Vk: previous potential
	% Ik: input current

	Km = 0.0041; 

	% km=0.14;
	Vth = 60000;  

	Vahp = 0;             

	Vk1=Vk+Ik-Vk*Km;

	if Vk1 >= Vth                                      
	Vk1 = Vahp;
	output=1;
	else 
	output=0;
	end
	*/
	m_value = 0.0; // for Population recording

	float Ke = 0.27;
	Ik = Ik -Ik*Ke;

	float Vk1 = Vk+Ik-Vk*Km;

	if(IsRecording())
	{
		if(Vk1>=Vth)
			m_storedValues.push_back(Vth);
		else
			m_storedValues.push_back(Vk1);
	}

	if(Vk1 >= Vth)
	{
		// spike

		Vk1 = Vahp;

		// send spike
		m_isNewEvent = true;
		float t = 0.0; // not used atm
		m_eventsOutgoing.push_back(this->CreateEvent(t));

		if(this->GetPopulation()->IsRecording() == true)
		{
			m_value = 1;
		}
	}

	Vk = Vk1;

	if(IsRecording())
	{
		m_storedValues.insert(m_storedValues.begin(),this->network()->GetCurrentTimeStep());
		m_recordedValues.push_back(m_storedValues);
		m_storedValues.clear();
	}
}


UnitModifier* UnitIF::CreateEvent(float time)
{
	// m_hypercolumnId ?
	UnitModifierSpike* e = new UnitModifierSpike(m_unitId,m_hypercolumnId,time);
	return e;
}


#if GSL_AVAILABLE==1

//#include "gsl_errno.h"
//#include "NetworkUnitsIF.h"

/* BeginDocumentation
Name: aeif_cond_alpha -  Conductance based exponential integrate-and-fire neuron model according to Brette and Gerstner (2005).

Description:
aeif_cond_alpha is the adaptive exponential integrate and fire neuron according to Brette and Gerstner (2005).
Synaptic conductances are modelled as alpha-functions.

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg solver with adaptive stepsize to integrate
the differential equation.

The membrane potential is given by the following differential equation:
C dV/dt= -g_L(V-E_L)+g_L*Delta_T*exp((V-V_T)/Delta_T)-g_e(t)(V-E_e) -g_i(t)(V-E_i)-w +I_e

and

tau_w * dw/dt= a(V-E_L) -W

Parameters: 
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m        double - Membrane potential in mV
  g_ex       double - Excitatory synaptic conductance in nS.
  dg_ex      double - First derivative of g_ex in nS/ms
  g_in       double - Inhibitory synaptic conductance in nS.
  dg_in      double - First derivative of g_in in nS/ms.
  w          double - Spike-adaptation current in pA.

Membrane Parameters:
  C_m        double - Capacity of the membrane in pF
  t_ref      double - Duration of refractory period in ms. 
  V_peak     double - Spike detection threshold in mV.
  V_reset    double - Reset value for V_m after a spike. In mV.
  E_L        double - Leak reversal potential in mV. 
  g_L        double - Leak conductance in nS.
  I_e        double - Constant external input current in pA.

Spike adaptation parameters:
  a          double - Subthreshold adaptation in nS.
  b          double - Spike-triggered adaptation in pA.
  Delta_T    double - Slope factor in mV
  tau_w      double - Adaptation time constant in ms
  V_t        double - Spike initiation threshold in mV (V_th can also be used for compatibility).

Synaptic parameters
  E_ex       double - Excitatory reversal potential in mV.
  tau_syn_ex double - Rise time of excitatory synaptic conductance in ms (alpha function).
  E_in       double - Inhibitory reversal potential in mV.
  tau_syn_in double - Rise time of the inhibitory synaptic conductance in ms (alpha function).

Author: Marc-Oliver Gewaltig

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, PotentialRequest, SynapticConductanceRequest, AeifWRequest

References: Brette R and Gerstner W (2005) Adaptive Exponential Integrate-and-Fire Model as 
            an Effective Description of Neuronal Activity. J
            Neurophysiol 94:3637-3642

SeeAlso: iaf_cond_alpha
*/
//
//

using namespace std;

extern "C"
 int aeif_cond_alpha_dynamics (double, const double y[], double f[], void* param)
 {
   // shorthand for class we work for
   //typedef nest::aeif_cond_alpha AEIF;
   
   // get parameters as reference  
   UnitadEIF::Parameters_* tmp =
     reinterpret_cast<UnitadEIF::Parameters_*>(param);
//   assert(tmp);
   UnitadEIF::Parameters_& p = *tmp;

   // shorthand for state variables
   const double& V     = std::min(0.0,y[UnitadEIF::V_M  ]); // fix
   const double& dg_ex = y[UnitadEIF::DG_EX];
   const double&  g_ex = y[UnitadEIF::G_EX ];
   const double& dg_in = y[UnitadEIF::DG_IN];
   const double&  g_in = y[UnitadEIF::G_IN ];
   const double& w     = y[UnitadEIF::W    ];

   const double I_syn_exc = g_ex * (V - p.E_ex);
   const double I_syn_inh = g_in * (V - p.E_in);
   const double I_spike = p.Delta_T * exp((V - p.V_th) / p.Delta_T);

   // dv/dt
   f[UnitadEIF::V_M  ] = ( -p.g_L *( (V-p.E_L) - I_spike) 
 	                   - I_syn_exc - I_syn_inh - w + p.I_e + p.I_stim) / p.C_m;

   f[UnitadEIF::DG_EX] = -dg_ex / p.tau_syn_ex;
   f[UnitadEIF::G_EX ] =  dg_ex - g_ex / p.tau_syn_ex; // Synaptic Conductance (nS)

   f[UnitadEIF::DG_IN] = -dg_in / p.tau_syn_in;
   f[UnitadEIF::G_IN ] =  dg_in - g_in / p.tau_syn_in; // Synaptic Conductance (nS)

   // Adaptation current w.
   f[UnitadEIF::W    ] = ( p.a * (V - p.E_L) - w ) / p.tau_w;

   return 0;//GSL_SUCCESS;
 }

UnitadEIF::Parameters_::Parameters_()
  : V_peak_    (  0.0    ),  // mV
    V_reset_   (-60.0    ),  // mV
    t_ref_     (  0.0    ),  // ms
    g_L        ( 30.0    ),  // nS
    C_m        (281.0    ),  // pF
    E_ex       (  0.0    ),  // mV
    E_in       (-85.0    ),  // mV
    E_L        (-70.6    ),  // mV
    Delta_T    (  2.0    ),  // mV
    tau_w      (144.0    ),  // ms
    a          (400.0    ),  // pS (nS?)
    b          ( 80.5    ),  // pA
    V_th       (-50.4    ),  // mV
    tau_syn_ex (  0.2    ),  // ms
    tau_syn_in (  2.0    ),  // ms
    I_e        (  0.0    )   // pA
{}


UnitadEIF::State_::State_(const Parameters_& p)
  : r_(0)
{
  y_[0] = p.E_L;
  for ( size_t i = 1 ; i < NSTATES ; ++i )
    y_[i] = 0;
}

UnitadEIF::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < NSTATES ; ++i )
    y_[i] = s.y_[i];
}

UnitadEIF::State_& UnitadEIF::State_::operator=(const State_& s)
{
//  assert(this != &s);  // would be bad logical error in program

  for ( size_t i = 0 ; i <NSTATES ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

UnitadEIF::UnitadEIF()
{
	m_unitType = "minicolumn"; // currently name for all comp units

/*	B_.spike_exc_.clear();       // includes resize
	B_.spike_inh_.clear();       // includes resize
	B_.currents_.clear();        // includes resize
	B_.potentials_.clear_data(); // includes resize
	B_.conductances_.clear_data(); // includes resize
	B_.aeif_ws_.clear_data();      // includes resize
	Archiving_Node::clear_history();
*/
/*	//B_.step_ = Time::get_resolution().get_ms();
	B_.step_ = this->network()->GetTimeResolution();//this->Population()->network()->GetSimulationResolution();

	// We must integrate this model with high-precision to obtain decent results
	B_.IntegrationStep_ = std::min(0.01, B_.step_);
*/
	static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;  // T

	if ( B_.s_ == 0 )
		B_.s_ = gsl_odeiv_step_alloc (T1, NSTATES);// T
	else 
		gsl_odeiv_step_reset(B_.s_);

	if ( B_.c_ == 0 )  
		B_.c_ = gsl_odeiv_control_yp_new (1e-6,1e-6);// T
	else
		gsl_odeiv_control_init(B_.c_, 1e-6, 1e-6, 0.0, 1.0);// T

	if ( B_.e_ == 0 )  
		B_.e_ = gsl_odeiv_evolve_alloc(NSTATES);
	else 
		gsl_odeiv_evolve_reset(B_.e_);

	B_.sys_.function  = aeif_cond_alpha_dynamics; 
	B_.sys_.jacobian  = NULL;
	B_.sys_.dimension =  NSTATES;
	B_.sys_.params    = reinterpret_cast<void*>(&P_);

	// change to correct values
	S_.y_[0] = P_.E_L; S_.y_[1] = S_.y_[2] =S_.y_[3] =S_.y_[4] =S_.y_[5] =0.0;
	S_.r_ = 0;
	V_.RefractoryCounts_ = 0; //?
	P_.I_stim = 0.0;

	V_.g0_ex_  = 1.0 * exp(1.0f) / P_.tau_syn_ex;
	V_.g0_in_  = 1.0 * exp(1.0f) / P_.tau_syn_in;
	//V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref_)).get_steps();
}

UnitadEIF::Buffers_::Buffers_()
  : s_(0),
    c_(0),
    e_(0)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

UnitadEIF::State_::State_()
{
}

UnitadEIF::~UnitadEIF()
{
//T
	if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
	if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
	if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

UnitModifier* UnitadEIF::CreateEvent(float time)
{
	UnitModifierSpike* e = new UnitModifierSpike(m_unitId,m_hypercolumnId,time);
	/*e->SetValue(1.0);
	e->SetTime(time);
	e->SetFromUnitId(m_unitId);*/
	//e->SetFromHypercolumnId(m_hypercolumnId);
	return e;
}

void UnitadEIF::update()
{
	double t = 0.0;

	// move to init, but then fail to change during runtime etc
	B_.step_ = this->network()->GetTimeResolution();//this->Population()->network()->GetSimulationResolution();
	B_.IntegrationStep_ = std::min(0.01, B_.step_);

	m_value = 0;

	float V0 = S_.y_[0];
	float V1 = S_.y_[1];
	float V2 = S_.y_[2];
	float V3 = S_.y_[3];
	float V4 = S_.y_[4];
	float V5 = S_.y_[5];

	if(S_.r_>0)
		--S_.r_;

	while ( t < B_.step_ )
	{
//T  
		V0 = S_.y_[0];
		V1 = S_.y_[1];
		V2 = S_.y_[2];
		V3 = S_.y_[3];
		V4 = S_.y_[4];
		V5 = S_.y_[5];
		
		const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
		  	   &B_.sys_,             // system of ODE
			   &t,                   // from t
			    B_.step_,            // to t <= step
			   &B_.IntegrationStep_, // integration step size
			    S_.y_);              // neuronal state


		if(B_.IntegrationStep_ < B_.step_/10.0)
		{
			bool b = false;
			B_.IntegrationStep_ = B_.step_/10.0;
		}
	  // spikes are handled inside the while-loop
	  // due to spike-driven adaptation

		if(S_.y_[V_M] != S_.y_[V_M] || S_.y_[W] != S_.y_[W] )
		{
			// currently handling errors with assuming spike (this is the cause) - change this so no error can occur
			//return;
			S_.y_[V_M] = P_.V_peak_;
			S_.y_[W] = V5;
		}

	  if ( S_.r_ > 0 )
	  {
			S_.y_[V_M] = P_.V_reset_;
	  }
	  else if ( S_.y_[V_M] >= P_.V_peak_ )
	  {
		  S_.y_[V_M]  = P_.V_reset_;
		  S_.y_[W]   += P_.b; // spike-driven adaptation
		  S_.r_       = 10;//V_.RefractoryCounts_;

		  if(this->Population()->IsRecording() == true)
		  {
			  m_value = 1;
		  }

		  m_isNewEvent = true;
		  m_eventsOutgoing.push_back(this->CreateEvent(t));
		  //return;
	  }

	  if(IsRecording())
	  {
		  if(S_.y_[V_M]>= P_.V_peak_)
			  m_storedValues.push_back(P_.V_peak_);
		  else
			  m_storedValues.push_back(S_.y_[V_M]);
	  }

	  //if ( status != GSL_SUCCESS )
      //  throw GSLSolverFailure(get_name(), status);
	}

	/*S_.y_[DG_EX] += B_.spike_exc_.get_value(lag) * V_.g0_ex_;
    S_.y_[DG_IN] += B_.spike_inh_.get_value(lag) * V_.g0_in_;
      
    // set new input current
    P_.I_stim = B_.currents_.get_value(lag);

    // voltage logging
    B_.potentials_.record_data(origin.get_steps()+lag, S_.y_[V_M]);
    B_.conductances_.record_data(origin.get_steps()+lag, 
				   std::pair<double_t, double_t>(S_.y_[G_EX], S_.y_[G_IN]));
    B_.aeif_ws_.record_data(origin.get_steps()+lag, S_.y_[W]);*/
}

void UnitadEIF::SimulateEventQueue()
{
	if(m_eventsIncoming.size()>0)
	{
		//vector<long> currentIncomingSpikes = m_eventsIncoming[0];
		vector<long> currentIncomingSpikes = m_eventsIncoming[0];
		
		if(currentIncomingSpikes.size()>0)
		{
			float totWeights = 0;

			m_eventsIncoming[0].clear();

			vector<Connection*> conns = this->Population()->GetIncomingConnections();
			vector<UnitModifier*> eus;
			vector<float> weights;

			// put connection vars in network instead of connection to avoid check
			for(int j=0;j<1;j++)//conns.size();j++)
			{
				for(int i=0;i<currentIncomingSpikes.size();i++)//m_eventsIncoming.size();i++)
				{
					UnitModifier* eu = network()->GetUnitModifierIncoming(currentIncomingSpikes[i]);

					//m_value += m_eventsIncoming[i]->GetValue() * conns[j]->GetWeight(m_eventsIncoming[i]->GetFromUnitId(),this->GetUnitId());
					bool exists = true;

					if(conns.size()>1) // if more than one incoming connection we need to know which connection it is coming from
					{
						exists = true;//network()->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());//conns[j]->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());
					}

					if(exists)
					{
						//totWeights += conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
						float weight = network()->GetWeight(eu->GetFromUnitId(),this->GetUnitId());//conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
						if(weight<0)
							S_.y_[DG_IN] += -weight * V_.g0_in_;
						else
							S_.y_[DG_EX] += weight * V_.g0_ex_;
						
						//weights.push_back(weight);
						//eus.push_back(eu);
						//eu->SetValue(eu->GetValue() * weight);
					}
				}
			}

			//if(S_.y_[DG_EX] > 300 * V_.g0_ex_)
			//	S_.y_[DG_EX] = 300 * V_.g0_ex_;
		}

		// elsewhere now
		//m_eventsIncoming.erase(m_eventsIncoming.begin(),m_eventsIncoming.begin()+1);
	}

	update();

	// let connections modify incoming
	int index=0;
	/*for(int i=0;i<m_preConnections.size();i++)
	{
		if(m_preConnections[i]->Pre()->IsLocal() == true)
		{
			Connection* c = m_preConnections[i];
			c->SimulateEvent(m_eventsIncoming[index]);
			index++;
		}
	}*/

	/*for(int i=0;i<m_eventsIncoming.size();i++)
	{
		Connection* c = m_hashIdConnection[m_eventsIncoming[i]->GetFromUnitId()];
		c->SimulateEvent(m_eventsIncoming[index]);
	}*/

	// weights



/*	vector<Connection*> conns = this->Population()->GetIncomingConnections();
	for(int j=0;j<conns.size();j++)
	for(int i=0;i<m_eventsIncoming.size();i++)
	{
		//m_value += m_eventsIncoming[i]->GetValue() * conns[j]->GetWeight(m_eventsIncoming[i]->GetFromUnitId(),this->GetUnitId());
		float weight = conns[j]->GetWeight(m_eventsIncoming[i]->GetFromUnitId(),this->GetUnitId());
		m_eventsIncoming[i]->SetValue(m_eventsIncoming[i]->GetValue() * weight);
	}

	m_eventsIncoming.clear();
	*/

	if(IsRecording())
	{
//		vector<float> values(1);
//		values[0] = m_storedValues;//S_.y_[V_M];
		m_storedValues.insert(m_storedValues.begin(),this->network()->GetCurrentTimeStep());
		m_recordedValues.push_back(m_storedValues);
		m_storedValues.clear();
	}
}


vector<vector<float> > UnitadEIF::GetValuesToRecord() 
{
	/*vector<vector<float> > f(1);
	vector<float> f2(1);
	f2[0] = S_.y_[V_M];//m_value;
	f[0] = f2;

	return f;*/
	return m_recordedValues;
}
/*
UnitIF::Parameters_::Parameters_()
  : C_      (250.0    ),  // pF
    Tau_    ( 10.0    ),  // ms
    tau_syn_(  2.0    ),  // ms
    TauR_   (  2.0    ),  // ms
    U0_     (-70.0    ),  // mV
    V_reset_(-70.0-U0_),  // mV, rel to U0_
    Theta_  (-55.0-U0_),  // mV, rel to U0_
    I_e_    (  0.0    )   // pA
{}

UnitIF::State_::State_()
  : y0_(0.0),
    y1_(0.0),
    y2_(0.0),
    y3_(0.0),  
    r_ (0)
{}
*/

/*void UnitIF::update(const long from, const long to)//Time const & origin, const long_t from, const long_t to)
{
//  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
//  assert(from < to);

  for ( long lag = from ; lag < to ; ++lag )
    {
      if ( S_.r_ == 0 )
	{
	  // neuron not refractory
	  S_.y3_ = V_.P30_*(S_.y0_ + P_.I_e_) + V_.P31_*S_.y1_ + V_.P32_*S_.y2_ + V_.P33_*S_.y3_;
	}
      else // neuron is absolute refractory
	--S_.r_;

      // alpha shape PSCs
      S_.y2_  = V_.P21_*S_.y1_ + V_.P22_ * S_.y2_;
      S_.y1_ *= V_.P11_;
    
      // Apply spikes delivered in this step: The spikes arriving at T+1 have an 
      // immediate effect on the state of the neuron
//      S_.y1_ += V_.PSCInitialValue_* B_.spikes_.get_value(lag);   
    
      // threshold crossing
      if (S_.y3_ >= P_.Theta_)
	{
	  S_.r_ = V_.RefractoryCounts_;
	  S_.y3_= P_.V_reset_; 
      
	  // A supra-threshold membrane potential should never be observable.
	  // The reset at the time of threshold crossing enables accurate integration
	  // independent of the computation step size, see [2,3] for details.   
	  
	  //set_spiketime(Time::step(origin.get_steps()+lag+1));
	  
	  //SpikeEvent se;
	  //network()->send(*this, se, lag);
	}

      // set new input current
//      S_.y0_ = B_.currents_.get_value(lag);

      // voltage logging
      //B_.potentials_.record_data(origin.get_steps()+lag, S_.y3_ + P_.U0_);
    }
}                          


UnitIF::Parameters_::Parameters_()
  : C_      (250.0    ),  // pF
    Tau_    ( 10.0    ),  // ms
    tau_syn_(  2.0    ),  // ms
    TauR_   (  2.0    ),  // ms
    U0_     (-70.0    ),  // mV
    V_reset_(-70.0-U0_),  // mV, rel to U0_
    Theta_  (-55.0-U0_),  // mV, rel to U0_
    I_e_    (  0.0    )   // pA
{}

UnitIF::State_::State_()
  : y0_(0.0),
    y1_(0.0),
    y2_(0.0),
    y3_(0.0),  
    r_ (0)
{}

UnitIF::UnitIF()
{
	m_unitType = "minicolumn";
	//const double h = Time::get_resolution().get_ms(); 
	double h = 0.01;

  // these P are independent
  V_.P11_ = V_.P22_ = exp(-h/P_.tau_syn_);
  V_.P33_ = exp(-h/P_.Tau_);
  V_.P21_ = h * V_.P11_;
  
  // these depend on the above. Please do not change the order.
  V_.P30_ = 1/P_.C_*(1-V_.P33_)*P_.Tau_;
  V_.P31_ = 1/P_.C_ * ((V_.P11_-V_.P33_)/(-1/P_.tau_syn_- -1/P_.Tau_)- h*V_.P11_)
    /(-1/P_.Tau_ - -1/P_.tau_syn_);
  V_.P32_ = 1/P_.C_*(V_.P33_-V_.P11_)/(-1/P_.Tau_ - -1/P_.tau_syn_);
  V_.PSCInitialValue_=1.0 * exp(1.0) /P_.tau_syn_;


  // TauR specifies the length of the absolute refractory period as 
  // a double_t in ms. The grid based iaf_neuron can only handle refractory
  // periods that are integer multiples of the computation step size (h).
  // To ensure consistency with the overall simulation scheme such conversion
  // should be carried out via objects of class nest::Time. The conversion 
  // requires 2 steps:
  //     1. A time object is constructed defining representation of 
  //        TauR in tics. This representation is then converted to computation time
  //        steps again by a strategy defined by class nest::Time.
  //     2. The refractory time in units of steps is read out get_steps(), a member
  //        function of class nest::Time.
  //
  // The definition of the refractory period of the iaf_neuron is consistent 
  // the one of iaf_psc_alpha_ps.
  //
  // Choosing a TauR that is not an integer multiple of the computation time 
  // step h will lead to accurate (up to the resolution h) and self-consistent
  // results. However, a neuron model capable of operating with real valued spike
  // time may exhibit a different effective refractory time.

  V_.RefractoryCounts_ = 1;//Time(Time::ms(P_.TauR_)).get_steps();
//  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

UnitModifier* UnitIF::CreateEvent()
{
	// should depend on the connection
	// this->m_pos
	UnitModifierGraded* e = new UnitModifierGraded;
	e->SetValue(1.0);
	e->SetFromUnitId(m_unitId);
	//e->SetFromHypercolumnId(m_hypercolumnId);

	return e;
}

void UnitIF::SimulateEventQueue()
{
	update(0,1);
}*/


#endif