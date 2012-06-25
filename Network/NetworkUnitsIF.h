#pragma once
#ifndef NETWORKUNITSIF_H
#define NETWORKUNITSIF_H


//#include "Network.h"
#include "NetworkUnits.h"

#if GSL_AVAILABLE==1

#include <gsl/gsl_odeiv.h>

#endif

class Unit;

using namespace std;


class UnitIF : public Unit
{
public:
	UnitIF();

	void SimulateEventQueue();
	UnitModifier* CreateEvent() { return CreateEvent(0.0); }
	UnitModifier* CreateEvent(float time);
	
	vector<vector<float> > GetValuesToRecord()
	{
		return m_recordedValues;
	}

	void SetUnitIdLocalHypercolumn(int id)
	{
		m_localIdHypercolumn = id;
	}

	void SetHypercolumnId(int id)
	{
		m_hypercolumnId = id;
	}

	void AddHypercolumn(Hypercolumn* h)
	{
		m_hypercolumns.push_back(h);
	}

private:

	vector<float> m_storedValues;
	int m_hypercolumnId; // parent hypercolumn id (used?)
	int m_localIdHypercolumn;
	vector<Hypercolumn*> m_hypercolumns;

	// naming from matlab script
	float Km,Vth,Vahp,Vk,Ik;
	
};
#if GSL_AVAILABLE==1

// Model brought in with small modifications from NEST
// Equation solving exactly as in NEST

class UnitadEIF : public Unit
{
public:
	UnitadEIF();
	~UnitadEIF();

	void SimulateEventQueue();
	void update(); // nest-retrieved fcn

	UnitModifier* CreateEvent() { return CreateEvent(0.0); }
	UnitModifier* CreateEvent(float time);

	enum Statevars
	{
		V_M   = 0,
		DG_EX    ,  // 1
		G_EX     ,  // 2
		DG_IN    ,  // 3
		G_IN     ,  // 4
		W        ,  // 5
		NSTATES
	};

	struct Parameters_ {
		double V_peak_;     //!< Spike detection threshold in mV
		double V_reset_;    //!< Reset Potential in mV
		double t_ref_;      //!< Refractory period in ms

		double g_L;         //!< Leak Conductance in nS
		double C_m;         //!< Membrane Capacitance in pF
		double E_ex;        //!< Excitatory reversal Potential in mV
		double E_in;        //!< Inhibitory reversal Potential in mV
		double E_L;         //!< Leak reversal Potential (aka resting potential) in mV
		double Delta_T;     //!< Slope faktor in ms.
		double tau_w;       //!< adaptation time-constant in ms.
		double a;           //!< Subthreshold adaptation in nS.
		double b;           //!< Spike-triggered adaptation in pA
		double V_th;        //!< Spike threshold in mV.
		//double t_ref;       //!< Refractory period in ms.
		double tau_syn_ex;  //!< Excitatory synaptic rise time.
		double tau_syn_in;  //!< Excitatory synaptic rise time.
		double I_e;         //!< Intrinsic current in pA.

		/** 
		* External input current from CurrentEvents.
		* This is not a parameter but a variable. It is still placed here, since
		* it needs to be passed to the iteration function. We thus avoid the need
		* of an additional wrapper structure. It is not revealed or manipulateable.
		*/
		double I_stim;      //!< External Stimulus in pA

		Parameters_();  //!< Sets default parameter values

		//void get(DictionaryDatum&) const;  //!< Store current values in dictionary
		//void set(const DictionaryDatum&);  //!< Set values from dicitonary
	};

	// functions for columnar structure - should merge this and minicolumn
	void SetUnitIdLocalHypercolumn(int id)
	{
		m_localIdHypercolumn = id;
	}

	void SetHypercolumnId(int id)
	{
		m_hypercolumnId = id;
	}
	void AddHypercolumn(Hypercolumn* h)
	{
		m_hypercolumns.push_back(h);
	}

	vector<vector<float> > GetValuesToRecord();

	void SetParameters(Parameters_ params)
	{
		P_=params;
	}

	void AddInput(float value) // special fcn, remove
	{
		S_.y_[DG_EX]+=value;
	}

private:

	vector<float> m_storedValues;
	int m_localIdHypercolumn;
	int m_hypercolumnId; // parent hypercolumn id
	vector<Hypercolumn*> m_hypercolumns;

	/**
	* Buffers of the model.
	*/
	struct Buffers_ {
		Buffers_(); //!<Sets buffer pointers to 0
		/** buffers and sums up incoming spikes/currents */
		//RingBuffer spike_exc_;
		//RingBuffer spike_inh_;
		//RingBuffer currents_;

		//AnalogDataLogger<PotentialRequest>           potentials_;
		//AnalogDataLogger<SynapticConductanceRequest> conductances_;
		//AnalogDataLogger<AeifWRequest>               aeif_ws_;

		/** GSL ODE stuff */
		gsl_odeiv_step*    s_;    //!< stepping function
		gsl_odeiv_control* c_;    //!< adaptive stepsize control function
		gsl_odeiv_evolve*  e_;    //!< evolution function
		gsl_odeiv_system   sys_;  //!< struct describing system

		// IntergrationStep_ should be reset with the neuron on ResetNetwork,
		// but remain unchanged during calibration. Since it is initialized with
		// step_, and the resolution cannot change after nodes have been created,
		// it is safe to place both here.
		double step_;           //!< step size in ms
		double   IntegrationStep_;//!< current integration time step, updated by GSL
	};

	    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {
		State_(); 
      double y_[NSTATES];  //!< neuron state, must be C-array for GSL solver
      int    r_;           //!< number of refractory steps remaining

      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      //void get(DictionaryDatum&) const;
      //void set(const DictionaryDatum&, const Parameters_&);
    };    

	/**
	* Internal variables of the model.
	*/
	struct Variables_ { 
		/** initial value to normalise excitatory synaptic conductance */
		double g0_ex_; 

		/** initial value to normalise inhibitory synaptic conductance */
		double g0_in_;    

		int    RefractoryCounts_;
	};

	Parameters_ P_; // should be replaced to a pointer
	Buffers_    B_;
	State_      S_;
	Variables_  V_;

	/*double V_peak_;     //!< Spike detection threshold in mV
	double V_reset_;    //!< Reset Potential in mV
	double t_ref_;      //!< Refractory period in ms

	double g_L;         //!< Leak Conductance in nS
	double C_m;         //!< Membrane Capacitance in pF
	double E_ex;        //!< Excitatory reversal Potential in mV
	double E_in;        //!< Inhibitory reversal Potential in mV
	double E_L;         //!< Leak reversal Potential (aka resting potential) in mV
	double Delta_T;     //!< Slope faktor in ms.
	double tau_w;       //!< adaptation time-constant in ms.
	double a;           //!< Subthreshold adaptation in nS.
	double b;           //!< Spike-triggered adaptation in pA
	double V_th;        //!< Spike threshold in mV.
	double t_ref;       //!< Refractory period in ms.
	double tau_syn_ex;  //!< Excitatory synaptic rise time.
	double tau_syn_in;  //!< Excitatory synaptic rise time.
	double I_e;         //!< Intrinsic current in pA.
	*/

	gsl_odeiv_step*    s_;    //!< stepping function
	gsl_odeiv_control* c_;    //!< adaptive stepsize control function
	gsl_odeiv_evolve*  e_;    //!< evolution function
	gsl_odeiv_system   sys_;  //!< struct describing system

	// IntergrationStep_ should be reset with the neuron on ResetNetwork,
	// but remain unchanged during calibration. Since it is initialized with
	// step_, and the resolution cannot change after nodes have been created,
	// it is safe to place both here.
	double step_;           //!< step size in ms
	double   IntegrationStep_;//!< current integration time step, updated by GSL

	//static const int NSTATES = 7; // ?

	/**
	* Enumeration identifying elements in state array State_::y_.
	* The state vector must be passed to GSL as a C array. This enum
	* identifies the elements of the vector. It must be public to be
	* accessible from the iteration function.

	enum Statevars
	{
	V_M   = 0,
	DG_EX    ,  // 1
	G_EX     ,  // 2
	DG_IN    ,  // 3
	G_IN     ,  // 4
	W        ,  // 5
	NSTATES
	};*/

	//double y_[NSTATES];  //!< neuron state, must be C-array for GSL solver

};


//public:
//
//	UnitIF()
//	{
//		m_unitType = "IF";
//		m_nrEventsSent = 0;
//		m_nrEventsReceived = 0;
//	}
//
//private:
//
//	void update(const long from, const long to);
//	float m_value;
//
///** 
//     * Independent parameters of the model. 
//     */
//    struct Parameters_ {
//      /** Membrane capacitance in pF. */
//      float C_;
//    
//      /** Membrane time constant in ms. */
//      float Tau_; 
//
//      /** Time constant of synaptic current in ms. */
//      float tau_syn_;
//      
//      /** Refractory period in ms. */
//      float TauR_;
//
//      /** Resting potential in mV. */
//      float U0_;
//
//      /** Reset value of the membrane potential, in mV.
//          @note Value is relative to resting potential U0_. */
//      float V_reset_;
//
//      /** Threshold in mV. 
//          @note Value is relative to resting potential U0_. */
//      float Theta_;
//
//      /** External current in pA */
//      float I_e_;
//
//      Parameters_();  //!< Sets default parameter values
//
//      //void get(DictionaryDatum&) const;  //!< Store current values in dictionary
//      //void set(const DictionaryDatum&);  //!< Set values from dicitonary
//    };
//
//    // ---------------------------------------------------------------- 
//
//    /**
//     * State variables of the model.
//     */
//    struct State_ {
//      float y0_; //!< Constant current
//      float y1_;  
//      float y2_;
//      float y3_; //!< This is the membrane potential RELATIVE TO RESTING POTENTIAL.
//
//      int    r_;  //!< number of refractory steps remaining
//
//      State_();  //!< Default initialization
//      
//      //void get(DictionaryDatum&, const Parameters_&) const;
//      //void set(const DictionaryDatum&, const Parameters_&);
//    };
//
//
//	    /**
//     * Buffers of the model.
//     */
//    struct Buffers_ {
//      /** buffers and summs up incoming spikes/currents */
//      RingBuffer spikes_;
//      RingBuffer currents_;
//
//      /**
//       * Buffer for membrane potential.
//       */
//    //  AnalogDataLogger<PotentialRequest> potentials_;
//    };
//
//	
//    /**
//     * Internal variables of the model.
//     */
//    struct Variables_ { 
//      /** Amplitude of the synaptic current.
//	        This value is chosen such that a post-synaptic potential with
//	        weight one has an amplitude of 1 mV.
//       */
//     float PSCInitialValue_;
//     int    RefractoryCounts_;  //!< refractory time in steps
//    
//     float P11_;   
//     float P21_;
//     float P22_;
//     float P31_;
//     float P32_;
//     float P30_;
//     float P33_;
//     
//   };
//
//	   /**
//    * @defgroup iaf_neuron_data
//    * Instances of private data structures for the different types
//    * of data pertaining to the model.
//    * @note The order of definitions is crucial: Moving Variables_
//    *       to the very end increases simulation time for brunel-2.sli
//    *       from 72s to 81s on a Mac, Intel Core 2 Duo 2.2GHz, g++ 4.0.1 -O3
//    * @{
//    */   
//   Parameters_ P_;
//   State_      S_;
//   Variables_  V_;
//   Buffers_    B_;
//   /** @} */
//};

#endif

#endif