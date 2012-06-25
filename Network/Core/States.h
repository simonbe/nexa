#pragma once

class States
{
public:

	States()
	{
		m_on = true;
	}

	bool IsOn()
	{
		return m_on;
	}

	void SwitchOnOff(bool on)
	{
		m_on = on;
	}

private:

	bool m_on;
};