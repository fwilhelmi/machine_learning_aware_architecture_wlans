# A Flexible Machine Learning-Aware Architecture for Future WLANs
### Authors
* [Francesc Wilhelmi](https://fwilhelmi.github.io/)
* [Sergio Barrachina-Muñoz](https://github.com/sergiobarra)
* [Boris Bellalta](http://www.dtic.upf.edu/~bbellalt/)
* [Cristina Cano](http://ccanobs.github.io/)
* [Anders Jonsson](http://www.tecn.upf.es/~jonsson/)
* [Vishnu Ram](https://www.researchgate.net/profile/Vishnu_Ov)

### Project description
Lots of hopes have been placed in Machine Learning (ML) as a key enabler of future wireless networks. By taking advantage of the large volumes of data generated by networks, ML is expected to deal with the ever-increasing complexity of networking problems. Unfortunately, current networking systems are not yet prepared for supporting the ensuing requirements of ML-based applications, especially for enabling procedures related to data collection, processing, and output distribution. This article points out the architectural requirements that are needed to pervasively include ML as part of future wireless networks operation. To this aim, we propose to adopt the International Telecommunications Union (ITU) unified architecture for 5G and beyond. Specifically, we look into Wireless Local Area Networks (WLANs), which, due to their nature, can be found in multiple forms, ranging from cloud-based to edge-computing-like deployments. Based on the ITU's architecture, we provide insights on the main requirements and the major challenges of introducing ML to the multiple modalities of WLANs.
### Repository description
This repository contains the LaTeX files used for the magazine article "A Flexible Machine Learning-Aware Architecture for Future WLANs", which has been sent to "IEEE Communications Magazine".

### Simulation Results
The results that are presented in this paper have been obtained from the IEEE 802.11 standard-compliant simulation framework available at [https://github.com/toniadame/WiFi_AP_Selection_Framework](https://github.com/toniadame/WiFi_AP_Selection_Framework). **Disclaimer:** the code provided in this repository has been integrated to the AP association framework.

In particular, we demonstrated the superiority of the ITU-T's architecture by showcasing the benefits of applying an ML-enabled AP (re)association solution. Specifically, we have used a Vanilla Neural Network (NN) approach, which has been trained with the data obtained from 100,000 random deployments. The simulation scenario comprises dense deployments with variable number of STAs and traffic load. The following Table summarizes the characteristics of the considered random deployments:

| **Parameter**    | **Value**        |
|------------------|------------------|
| Number of STAs   | 1-24             |
| Load of STAs     | 1-15 Mbps        |
| Type of traffic  | Uplink           |
| Location of STAs | Randomly uniform |
| Location of APs  | Fixed            |

An example of a random deployment is illustrated in the following Figure:

<img src="https://github.com/fwilhelmi/machine_learning_aware_architecture_wlans/blob/master/Other%20resources/use_case_ap_selection/random_sta_deployment.png" alt="Random deployment"
	title="Random deployment" width="500" />

The simulation parameters are as follows:

|     | **Parameter**               | **Value**                                                                                                                         |
|-----|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| PHY | Central frequency           | 5 GHz                                                                                                                             |
|     | Path-loss model             | See [IEEE 802.11ax Residential scenario](https://mentor.ieee.org/802.11/dcn/14/11-14-0980-16-00ax-simulation-scenarios.docx) |
|     | Tx gain / Rx gain           | 0 / O dB                                                                                                                          |
|     | Background noise            | -95 dBm                                                                                                                           |
|     | Legacy OFDM symbol duration | 4 micro s                                                                                                                          |
|     | OFDM symbol duration        | 16 micro s                                                                                                                       |
|     | No. of spatial streams      | 2                                                                                                                                 |
|     | Tx power                    | 20 dBm                                                                                                                            |
|     | Rx sensitivity              | -90 dBm                                                                                                                           |
| MAC | SIFS / DIFS duration        | 16 / 34 micro s                                                                                                                 |
|     | RTS / CTS length            | 160 / 112 bits                                                                                                                    |
|     | MH / SF / TB / MD length    | 320 / 16 / 18 / 32 bits                                                                                                           |
|     | Allowed data rates          | 14.4 to 173.4 Mbps                                                                                                                |

To generate the dataset (see [output_stas.csv](https://github.com/fwilhelmi/machine_learning_aware_architecture_wlans/blob/master/Other%20resources/use_case_ap_selection/output_stas.csv)), we have employed function [generate_data_set.m](https://github.com/fwilhelmi/machine_learning_aware_architecture_wlans/blob/master/Other%20resources/use_case_ap_selection/generate_data_set.m). The training procedure that generates our throughput prediction function is found at [neural_net_train.m](https://github.com/fwilhelmi/machine_learning_aware_architecture_wlans/blob/master/Other%20resources/use_case_ap_selection/neural_net_train.m) (the output used in the paper is [nn_function.mat](https://github.com/fwilhelmi/machine_learning_aware_architecture_wlans/blob/master/Other%20resources/use_case_ap_selection/nn_function.mat). The features considered for training are:
1. Rate STA [bps]: the rate at which a given STA can transmit to a specific AP (based on the maximum allowed MCS).
2. Load STA [bps]: the traffic load generated at the STA, based on different application requirements.
3. Delivery ratio AP [%]: percentage of load that a given AP is able to process before a new (re)association.
4. Load AP [bps]: amount of load that the AP is serving before a new (re)association.

Our results compare the standard Strongest Signal First (SSF) with the NN approach. In total, 50 random deployments are considered for averaging purposes. The following Figure shows the obtained results for a variable number of STAs:

<img src="https://github.com/fwilhelmi/machine_learning_aware_architecture_wlans/blob/master/Other%20resources/use_case_ap_selection/results_use_case.png" alt="Results SSF vs NN"
	title="Results SSF vs NN" width="500" />

As shown, the NN approach improves the average satisfaction obtained in all the deployments. Apart from that, it provides stability and increases fairness.

### Contribute

If you want to contribute, please contact to [francisco.wilhelmi@upf.edu](francisco.wilhelmi@upf.edu)
