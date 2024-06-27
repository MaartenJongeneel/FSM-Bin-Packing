#include "FSMSequenceOne_Initial.h"

#include "../FSMSequenceOne.h"

void FSMSequenceOne_Initial::configure(const mc_rtc::Configuration & config) {}

void FSMSequenceOne_Initial::start(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<FSMSequenceOne &>(ctl_);
}

bool FSMSequenceOne_Initial::run(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<FSMSequenceOne &>(ctl_);
  output("OK");
  return true;
}

void FSMSequenceOne_Initial::teardown(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<FSMSequenceOne &>(ctl_);
}

EXPORT_SINGLE_STATE("FSMSequenceOne_Initial", FSMSequenceOne_Initial)
