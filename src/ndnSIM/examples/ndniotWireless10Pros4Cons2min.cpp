/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/**
 * Copyright (c) 2011-2015  Regents of the University of California.
 *
 * This file is part of ndnSIM. See AUTHORS for complete list of ndnSIM authors and
 * contributors.
 *
 * ndnSIM is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * ndnSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ndnSIM, e.g., in COPYING.md file.  If not, see <http://www.gnu.org/licenses/>.
 **/

// ndn-test.cpp
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
//#include "ns3/internet-module.h"

#include "ns3/ndnSIM-module.h"


//#include "ns3/core-module.h"
//#include "ns3/network-module.h"
//#include "ns3/point-to-point-module.h"
//#include "ns3/ndnSIM-module.h"

#include <sys/time.h>
#include "ns3/ndnSIM/utils/mem-usage.hpp"
#include "ns3/ndnSIM/model/cs/ndn-content-store.hpp"
#include "ns3/ndnSIM/utils/mem-usage.hpp"

namespace ns3 {

/**
 * This scenario simulates a very simple network topology:
 *
 *
 *      +----------+     10000Mbps   +----------+
 *      | consumer |  <------------> | producer |
 *      +----------+       10ms      +----------+
 *
 *
 *     NS_LOG=ndn.Consumer ./waf --run ndn-test
 */

class Tester {
public:
  Tester()
    : m_csSize(300)
    , m_payload("1024")
    , m_interestRate(500)
    , m_shouldEvaluatePit(false)
    , m_simulationTime(Seconds(60000) / m_interestRate)
  {
  }

  int
  run(int argc, char* argv[]);

  void
  printHeader(std::ostream& os);

  void
  printStats(std::ostream& os, Time nextPrintTime, double beginRealTime);

private:
  std::string m_oldContentStore;
  size_t m_csSize;
  std::string m_payload;
  double m_interestRate;
  bool m_shouldEvaluatePit;
  std::string m_strategy;
  double m_initialOverhead;
  double m_distance;
  Time m_simulationTime;
};

void
Tester::printHeader(std::ostream& os)
{
  m_initialOverhead = MemUsage::Get() / 1024.0 / 1024.0;
  os << "SimulationTime"
     << "\t"
     << "RealTime"
     << "\t"
     << "NumberOfInterests (total)"
     << "\t"
     << "NumberOfInterests (per real time)"
     << "\n";
}

void
Tester::printStats(std::ostream& os, Time nextPrintTime, double beginRealTime)
{
  ::timeval t;
  gettimeofday(&t, NULL);
  double realTime = t.tv_sec + (0.000001 * (unsigned)t.tv_usec) - beginRealTime;
  Time simTime = Simulator::Now();

  os << simTime << "\t";
  os << realTime << "\t";

  double interestCount = m_interestRate * simTime.ToDouble(Time::S);
  double nInterestsPerSec = interestCount / realTime;

  os << interestCount << "\t" << nInterestsPerSec << "\t";

  uint64_t pitCount = 0;
  uint64_t nameTreeCount = 0;
  uint64_t csCount = 0;
  uint64_t droppedPacketCounts = 0;
  for (NodeList::Iterator node = NodeList::Begin(); node != NodeList::End(); node++) {

    auto pitSize = (*node)->GetObject<ndn::L3Protocol>()->getForwarder()->getPit().size();
    if (pitSize != 0)
      pitCount += pitSize;

    auto droppedPacketSize = (*node)->GetObject<ndn::L3Protocol>()->getForwarder()->getPit().missedInt();
    if (droppedPacketSize != 0)
      droppedPacketCounts += droppedPacketSize;

    auto nameTreeSize = (*node)->GetObject<ndn::L3Protocol>()->getForwarder()->getNameTree().size();
    if (nameTreeSize != 0)
      nameTreeCount += nameTreeSize;

    

    if (true != true) {
      Ptr<ndn::ContentStore> cs = (*node)->GetObject<ndn::ContentStore>();
      if (cs != 0)
        csCount += cs->GetSize();
    }
    else {
      auto csSize = (*node)->GetObject<ndn::L3Protocol>()->getForwarder()->getCs().size();
      if (csSize != 0)
        csCount += csSize;
    }
  }

  os << "pit" << ": " << pitCount << "\t";
  os << "cs" << ": " << csCount << "\t";
  os << "nametree" << ": " << nameTreeCount << "\t";
  os << "Drop packet: " << droppedPacketCounts << "\t";

  os << "Mem: " << MemUsage::Get() / 1024.0 / 1024.0 << "MiB\n";

  if ((simTime + nextPrintTime) >= m_simulationTime) {
    double finalOverhead = MemUsage::Get() / 1024.0 / 1024.0;
    if (m_shouldEvaluatePit) {
      if (pitCount != 0) {
        os << "Approximate memory overhead per PIT entry:"
           <<  1000 * (finalOverhead - m_initialOverhead) / pitCount << "KiB\n";
      }
      else {
        os << "`The number of PIT entries is equal to zero\n";
      }
    }
    else {
      if (csCount != 0) {
        os << "Approximate memory overhead per CS entry:"
           <<  1000 * (finalOverhead - m_initialOverhead) / csCount << "KiB\n";
      }
      else {
        os << "The number of CS entries is equal to zero\n";
      }
    }
  }

  Simulator::Schedule(nextPrintTime, &Tester::printStats, this, ref(os), nextPrintTime,
                      beginRealTime);
}

int
Tester::run(int argc, char* argv[])
{
  // setting default parameters for PointToPoint links and channels
  std::string phyMode ("DsssRate11Mbps");
  double rss = -80;  // -dBm
  uint32_t packetSize = 1000; // bytes
  uint32_t numPackets = 1;
  double interval = 1.0; // seconds
  bool verbose = false;

  CommandLine cmd;

  cmd.AddValue ("phyMode", "Wifi Phy mode", phyMode);
  cmd.AddValue ("rss", "received signal strength", rss);
  cmd.AddValue ("packetSize", "size of application packet sent", packetSize);
  cmd.AddValue ("numPackets", "number of packets generated", numPackets);
  cmd.AddValue ("interval", "interval (seconds) between packets", interval);
  cmd.AddValue ("verbose", "turn on all WifiNetDevice log components", verbose);


  cmd.AddValue("old-cs", "Old content store to use "
                         "(e.g., ns3::ndn::cs::Lru, ns3::ndn::cs::Lfu, ...)",
               m_oldContentStore);
  cmd.AddValue("cs-size", "Maximum number of cached packets per node", m_csSize);
  cmd.AddValue("rate", "Interest rate", m_interestRate);
  cmd.AddValue("pit", "Perform PIT evaluation if this parameter is true",
               m_shouldEvaluatePit);
  cmd.AddValue("strategy", "Choose forwarding strategy "
                           "(e.g., /localhost/nfd/strategy/multicast, "
                           "/localhost/nfd/strategy/best-route, ...) ",
               m_strategy);
  cmd.AddValue("sim-time", "Simulation time", m_simulationTime);
  cmd.Parse(argc, argv);

  cmd.Parse (argc, argv);
  // Convert to time object
  Time interPacketInterval = Seconds (interval);

  // disable fragmentation for frames below 2200 bytes
  Config::SetDefault ("ns3::WifiRemoteStationManager::FragmentationThreshold", StringValue ("2200"));
  // turn off RTS/CTS for frames below 2200 bytes
  Config::SetDefault ("ns3::WifiRemoteStationManager::RtsCtsThreshold", StringValue ("2200"));
  // Fix non-unicast data rate to be the same as that of unicast
  Config::SetDefault ("ns3::WifiRemoteStationManager::NonUnicastMode", 
                      StringValue (phyMode));

  NodeContainer c;
  c.Create (17);
  
  std::string prefixP1 = "/prefix/texst/gff/hg/g1";
  std::string prefixP2 = "/prefix/texst/gff/hg/g2";
  std::string prefixP3 = "/prefix/texst/gff/hg/g3"; 
  std::string prefixP4 = "/prefix/texst/gff/hg/g4";
  std::string prefixP5 = "/prefix/texst/gff/hg/g5";
  std::string prefixP6 = "/prefix/texst/gff/hg/g6";
  std::string prefixP7 = "/prefix/texst/gff/hg/g7";
  std::string prefixP8 = "/prefix/texst/gff/hg/g8";
  std::string prefixP9 = "/prefix/texst/gff/hg/g9";
  std::string prefixP10 = "/prefix/texst/gff/hg/g10";

  WifiHelper wifi;
  if (verbose)
    {
      wifi.EnableLogComponents ();  // Turn on all Wifi logging
    }
  wifi.SetStandard (WIFI_PHY_STANDARD_80211n_2_4GHZ);

  YansWifiPhyHelper wifiPhy =  YansWifiPhyHelper::Default ();
  // This is one parameter that matters when using FixedRssLossModel
  // set it to zero; otherwise, gain will be added
  wifiPhy.Set ("RxGain", DoubleValue (0) ); 
  // ns-3 supports RadioTap and Prism tracing extensions for 802.11b
  wifiPhy.SetPcapDataLinkType (YansWifiPhyHelper::DLT_IEEE802_11_RADIO); 

  YansWifiChannelHelper wifiChannel = YansWifiChannelHelper::Default ();
  //YansWifiChannelHelper wifiChannel;
  wifiChannel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
  wifiChannel.AddPropagationLoss ("ns3::NakagamiPropagationLossModel");
 // wifiChannel.AddPropagationLoss ("ns3::RangePropagationLossModel", doubleValue(50));
  // The below FixedRssLossModel will cause the rss to be fixed regardless
  // of the distance between the two stations, and the transmit power
 // wifiChannel.AddPropagationLoss ("ns3::FixedRssLossModel","Rss",DoubleValue (rss));
 // wifiChannel.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
  wifiPhy.SetChannel (wifiChannel.Create ());

  // Add a non-QoS upper mac, and disable rate control
  NqosWifiMacHelper wifiMac = NqosWifiMacHelper::Default ();
  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
                                "DataMode",StringValue (phyMode),
                                "ControlMode",StringValue (phyMode));
  // Set it to adhoc mode
  wifiMac.SetType ("ns3::AdhocWifiMac");
  NetDeviceContainer devices = wifi.Install (wifiPhy, wifiMac, c);

  // Note that with FixedRssLossModel, the positions below are not 
  // used for received signal strength. 
  MobilityHelper mobility;
  Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator> ();
  positionAlloc->Add (Vector (-120.0, 150.0, 0.0));//0
  positionAlloc->Add (Vector (-120.0, 80.0, 0.0));//1
  positionAlloc->Add (Vector (-120.0, 10.0, 0.0));//2
  positionAlloc->Add (Vector (-190.0, 80.0, 0.0));//3
  positionAlloc->Add (Vector (-20.0, 80.0, 0.0));//4
  positionAlloc->Add (Vector (70.0, 80.0, 0.0));//5
  positionAlloc->Add (Vector (150.0, 80.0, 0.0));//6
  positionAlloc->Add (Vector (210.0, 10.0, 0.0));//7
  positionAlloc->Add (Vector (210.0, 150.0, 0.0));//8
  positionAlloc->Add (Vector (280.0, 150.0, 0.0));//9
  positionAlloc->Add (Vector (280.0, 80.0, 0.0));//10
  positionAlloc->Add (Vector (280.0, 10.0, 0.0));//11
  positionAlloc->Add (Vector (350.0, 150.0, 0.0));//12
  positionAlloc->Add (Vector (350.0, 80.0, 0.0));//13
  positionAlloc->Add (Vector (350.0, 10.0, 0.0));//14
  positionAlloc->Add (Vector (420.0, 150.0, 0.0));//15
  positionAlloc->Add (Vector (420.0, 10.0, 0.0));//16
  mobility.SetPositionAllocator (positionAlloc);
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobility.Install (c);

  ////////////////
  // 1. Install Wifi
  NetDeviceContainer wifiNetDevices = wifi.Install(wifiPhy, wifiMac, c);

  // 2. Install Mobility model
  mobility.Install(c);

  // 3. Install NDN stack
  //NS_LOG_INFO("Installing NDN stack");
  ndn::StackHelper ndnHelper;
  // ndnHelper.AddNetDeviceFaceCreateCallback (WifiNetDevice::GetTypeId (), MakeCallback
  // (MyNetDeviceFaceCallback));
   if (!m_oldContentStore.empty()) {
    ndnHelper.SetOldContentStore(m_oldContentStore, "MaxSize", std::to_string(m_csSize));
   }
  //ndnHelper.SetDefaultRoutes(true);
  ndnHelper.Install(c);

  // Set BestRoute strategy
  ndn::StrategyChoiceHelper::Install(c, "/", "/localhost/nfd/strategy/best-route");

  // 4. Set up applications
  //NS_LOG_INFO("Installing Applications");
  ndn::GlobalRoutingHelper ndnGlobalRoutingHelper;
  ndnGlobalRoutingHelper.InstallAll();

  ndn::AppHelper consumerHelper("ns3::ndn::ConsumerCbr");
  consumerHelper.SetPrefix(prefixP1);
  consumerHelper.SetPuid("A1");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(0));

 // ndn::AppHelper consumerHelper("ns3::ndn::ConsumerCbr");
  consumerHelper.SetPrefix(prefixP2);
  consumerHelper.SetPuid("A2");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(0));

  consumerHelper.SetPrefix(prefixP3);
  consumerHelper.SetPuid("A3");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(0));
  
  consumerHelper.SetPrefix(prefixP4);
  consumerHelper.SetPuid("A4");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(1));

  consumerHelper.SetPrefix(prefixP5);
  consumerHelper.SetPuid("A5");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(1));   
 
  consumerHelper.SetPrefix(prefixP6);
  consumerHelper.SetPuid("A6");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(1));

  consumerHelper.SetPrefix(prefixP7);
  consumerHelper.SetPuid("A7");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(2));   

  consumerHelper.SetPrefix(prefixP8);
  consumerHelper.SetPuid("A8");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(2));      

  consumerHelper.SetPrefix(prefixP9);
  consumerHelper.SetPuid("A9");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(3));   

  consumerHelper.SetPrefix(prefixP10);
  consumerHelper.SetPuid("A10");
  consumerHelper.SetAttribute("Frequency", DoubleValue(m_interestRate));
  consumerHelper.Install(c.Get(3));


  ndn::AppHelper producerHelper("ns3::ndn::Producer");
  producerHelper.SetPrefix(prefixP1);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(7));

  producerHelper.SetPrefix(prefixP2);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(8));

  producerHelper.SetPrefix(prefixP3);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(9));
   
  producerHelper.SetPrefix(prefixP4);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(10));

  producerHelper.SetPrefix(prefixP5);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(11));

  producerHelper.SetPrefix(prefixP6);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(12));

  producerHelper.SetPrefix(prefixP7);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(13));

  producerHelper.SetPrefix(prefixP8);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(14));

  producerHelper.SetPrefix(prefixP9);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(15));
  
  producerHelper.SetPrefix(prefixP10);
  producerHelper.SetAttribute("PayloadSize", StringValue(m_payload));
  producerHelper.Install(c.Get(16));

  ndnGlobalRoutingHelper.AddOrigins(prefixP1, c.Get(7));
  ndnGlobalRoutingHelper.AddOrigins(prefixP2, c.Get(8));
  ndnGlobalRoutingHelper.AddOrigins(prefixP3, c.Get(9));
  ndnGlobalRoutingHelper.AddOrigins(prefixP4, c.Get(10));
  ndnGlobalRoutingHelper.AddOrigins(prefixP5, c.Get(11));
  ndnGlobalRoutingHelper.AddOrigins(prefixP6, c.Get(12));
  ndnGlobalRoutingHelper.AddOrigins(prefixP7, c.Get(13));
  ndnGlobalRoutingHelper.AddOrigins(prefixP8, c.Get(14));
  ndnGlobalRoutingHelper.AddOrigins(prefixP9, c.Get(15));
  ndnGlobalRoutingHelper.AddOrigins(prefixP10, c.Get(16));
  
  //ndnGlobalRoutingHelper.AddOrigins("/test/prefix", c.Get(0));

  ndn::GlobalRoutingHelper::CalculateRoutes();
 /* ndnGlobalRoutingHelper.AddOrigins("/test/prefix", producer3);
  ndnGlobalRoutingHelper.AddOrigins(prefixP4, producer4);
  ndnGlobalRoutingHelper.AddOrigins(prefixP5, producer5);*/



  Simulator::Stop(m_simulationTime);

  struct ::timeval t;
  gettimeofday(&t, NULL);
  double beginRealTime = t.tv_sec + (0.000001 * (unsigned)t.tv_usec);
  Simulator::Schedule(Seconds(0), &Tester::printHeader, this, ref(std::cout));
  Simulator::Schedule(m_simulationTime / 200, &Tester::printStats, this, ref(std::cout),
                      m_simulationTime / 200, beginRealTime);

  L2RateTracer::InstallAll("drop-trace2.txt", Seconds(0.5));
  Simulator::Run();
  Simulator::Destroy();

  return 0;
}

} // namespace ns3

int
main(int argc, char* argv[])
{
  ns3::Tester tester;
  return tester.run(argc, argv);
}
