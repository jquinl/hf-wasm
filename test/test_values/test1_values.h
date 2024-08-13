
#define cgto_number 6
#define cgto_n 6
#define test_at_num 2

const int test_types[2] = {1,8};
const int test_pos[2][3] ={ {1.0,1.0,1.0},{9.0,9.0,9.0}};

 float S_test[cgto_number * cgto_number] = { 
    0.9999999908897999, 8.485806523168149e-16, 4.8976976063215463e-11, -8.326316139686609e-11, -8.326316139686609e-11, -8.326316139686609e-11, 
8.485806523168148e-16, 1.0000000165470366, 0.2367039416499711, 6.311319040196644e-16, 6.311319040196644e-16, 6.311319040196644e-16, 
4.8976976063215463e-11, 0.23670394164997102, 1.0000000268753377, 0.0, 0.0, 0.0, 
-8.32631613968661e-11, 6.311319040196645e-16, 0.0, 1.0000000195212955, 0.0, 0.0, 
-8.32631613968661e-11, 6.311319040196645e-16, 0.0, 0.0, 1.0000000195212955, 0.0, 
-8.32631613968661e-11, 6.311319040196645e-16, 0.0, 0.0, 0.0, 1.0000000195212957};

const float T_test[cgto_number * cgto_number] = { 0.7600318766425666, -8.47554861947802e-15, -2.4016439128066893e-10, 3.8904416156913516e-10, 3.8904416156913516e-10, 3.8904416156913516e-10, 
-8.475548619478104e-15, 29.003200425456576, -0.16801094296420926, 3.1381112690203975e-15, 3.1381112690203975e-15, 3.1381112690203975e-15, 
-2.4016439128066903e-10, -0.16801094296421026, 0.8081279766490592, 0.0, 0.0, 0.0, 
3.8904416156913526e-10, -3.833078014107721e-15, 0.0, 2.528731247558872, 0.0, 0.0, 
3.8904416156913526e-10, -3.833078014107721e-15, 0.0, 0.0, 2.5287312475588717, 0.0, 
3.890441615691352e-10, -3.833078014107721e-15, 0.0, 0.0, 0.0, 2.5287312475588726};

const float V_test[cgto_number * cgto_number] = { -1.8039639858790375, -2.030645113498314e-14, -9.720447416483193e-11, 1.5774705470090662e-10, 1.5774705470090662e-10, 1.5774705470090662e-10, 
-2.030645113498314e-14, -60.690624693321304, -7.200162236201033, 0.00015273299347455679, 0.00015273299347455679, 0.00015273299347455679, 
-9.720447416483193e-11, -7.200162236201034, -9.124087702333487, 0.0019280277232065458, 0.0019280277232065458, 0.0019280277232065458, 
1.5774705470090662e-10, 0.00015273299347455676, 0.001928027723206546, -9.05162368694427, -0.00011167016945184659, -0.00011167016945184659, 
1.5774705470090662e-10, 0.00015273299347455676, 0.001928027723206546, -0.00011167016945184659, -9.05162368694427, -0.00011167016945184659, 
1.5774705470090662e-10, 0.00015273299347455676, 0.001928027723206546, -0.00011167016945184659, -0.00011167016945184659, -9.05162368694427};

const float EE_test[cgto_number * cgto_number * cgto_number * cgto_number] = { 0.7746059298062673, 6.259055533316608e-17, 5.102602193945161e-12, -8.86604034846065e-12, -8.86604034846065e-12, -8.86604034846065e-12, 6.259055533316608e-17, 0.07216878418541066, 0.017082635398105808, -0.00015273299209336394, -0.00015273299209336394, -0.00015273299209336394, 5.10260219394516e-12, 0.017082635398105804, 0.07216878493079165, -0.0019280277056418272, -0.0019280277056418272, -0.0019280277056418272, -8.866040348460652e-12, -0.0001527329920933639, -0.0019280277056418274, 0.07216878440005935, 0.00011167016843450898, 0.00011167016843450898, -8.866040348460652e-12, -0.0001527329920933639, -0.0019280277056418274, 0.00011167016843450898, 0.07216878440005935, 0.00011167016843450898, -8.866040348460652e-12, -0.0001527329920933639, -0.0019280277056418274, 0.00011167016843450898, 0.00011167016843450898, 0.07216878440005935, 
6.25905553331661e-17, 1.6040631582876427e-30, 1.0513263897152411e-26, -1.6929642664672614e-26, -1.6929642664672614e-26, -1.6929642664672614e-26, 1.6040631582876427e-30, 2.1844791426148358e-15, 4.25537987485138e-16, -3.881818184146135e-17, -3.881818184146135e-17, -3.881818184146135e-17, 1.0513263897152411e-26, 4.2553798748513794e-16, 8.920180666433653e-16, -8.92478658317932e-17, -8.92478658317932e-17, -8.92478658317932e-17, -1.6929642664672612e-26, -3.881818184146136e-17, -8.92478658317932e-17, 8.922618031428675e-16, 1.669299767318975e-17, 1.669299767318975e-17, -1.6929642664672612e-26, -3.881818184146136e-17, -8.92478658317932e-17, 1.669299767318975e-17, 8.922618031428675e-16, 1.669299767318975e-17, -1.6929642664672612e-26, -3.881818184146136e-17, -8.92478658317932e-17, 1.669299767318975e-17, 1.669299767318975e-17, 8.922618031428675e-16, 
5.1026021939451584e-12, 1.051326389715241e-26, 1.417481549327364e-21, -2.4091087190047952e-21, -2.4091087190047952e-21, -2.4091087190047952e-21, 1.051326389715241e-26, 1.1512704837563123e-11, 2.7250778192552604e-12, -7.946365155061487e-14, -7.946365155061487e-14, -7.946365155061487e-14, 1.417481549327364e-21, 2.7250778192552604e-12, 1.1507004690826273e-11, -9.980982538000196e-13, -9.980982538000196e-13, -9.980982538000196e-13, -2.4091087190047952e-21, -7.946365155061488e-14, -9.980982538000196e-13, 1.1506575748930733e-11, 1.8554495532496332e-13, 1.8554495532496332e-13, -2.4091087190047952e-21, -7.946365155061488e-14, -9.980982538000196e-13, 1.8554495532496332e-13, 1.1506575748930733e-11, 1.8554495532496332e-13, -2.4091087190047952e-21, -7.946365155061488e-14, -9.980982538000196e-13, 1.8554495532496332e-13, 1.8554495532496332e-13, 1.1506575748930733e-11, 
-8.866040348460655e-12, -1.6929642664672612e-26, -2.4091087190047956e-21, 4.196583276467744e-21, 4.0939607446465595e-21, 4.0939607446465595e-21, -1.692964266467261e-26, -1.8610120567131596e-11, -4.405077201779598e-12, 1.3582477533571906e-13, 1.150745348691181e-13, 1.150745348691181e-13, -2.4091087190047956e-21, -4.4050772017795975e-12, -1.8605515587788372e-11, 1.7083991283589316e-12, 1.4491112343398354e-12, 1.4491112343398354e-12, 4.196583276467743e-21, 1.3582477533571906e-13, 1.7083991283589318e-12, -1.866996636040683e-11, -2.891914359470458e-13, -2.891914359470458e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.891914359470458e-13, -1.8572666436284773e-11, -2.405414738860167e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.891914359470458e-13, -2.405414738860167e-13, -1.8572666436284773e-11, 
-8.866040348460655e-12, -1.6929642664672612e-26, -2.4091087190047956e-21, 4.0939607446465595e-21, 4.196583276467744e-21, 4.0939607446465595e-21, -1.692964266467261e-26, -1.8610120567131596e-11, -4.405077201779598e-12, 1.150745348691181e-13, 1.3582477533571906e-13, 1.150745348691181e-13, -2.4091087190047956e-21, -4.4050772017795975e-12, -1.8605515587788372e-11, 1.4491112343398354e-12, 1.7083991283589316e-12, 1.4491112343398354e-12, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -1.857266643628477e-11, -2.891914359470458e-13, -2.405414738860167e-13, 4.196583276467743e-21, 1.3582477533571906e-13, 1.7083991283589318e-12, -2.891914359470458e-13, -1.866996636040683e-11, -2.891914359470458e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.405414738860167e-13, -2.891914359470458e-13, -1.8572666436284773e-11, 
-8.866040348460655e-12, -1.6929642664672612e-26, -2.4091087190047956e-21, 4.0939607446465595e-21, 4.0939607446465595e-21, 4.196583276467744e-21, -1.692964266467261e-26, -1.8610120567131596e-11, -4.405077201779598e-12, 1.150745348691181e-13, 1.150745348691181e-13, 1.3582477533571906e-13, -2.4091087190047956e-21, -4.4050772017795975e-12, -1.8605515587788372e-11, 1.4491112343398354e-12, 1.4491112343398354e-12, 1.7083991283589316e-12, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -1.857266643628477e-11, -2.405414738860167e-13, -2.891914359470458e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.405414738860167e-13, -1.857266643628477e-11, -2.891914359470458e-13, 4.196583276467743e-21, 1.3582477533571906e-13, 1.7083991283589318e-12, -2.891914359470458e-13, -2.891914359470458e-13, -1.866996636040683e-11, 
6.25905553331661e-17, 1.6040631582876427e-30, 1.0513263897152411e-26, -1.6929642664672614e-26, -1.6929642664672614e-26, -1.6929642664672614e-26, 1.6040631582876427e-30, 2.1844791426148358e-15, 4.25537987485138e-16, -3.881818184146135e-17, -3.881818184146135e-17, -3.881818184146135e-17, 1.0513263897152411e-26, 4.2553798748513794e-16, 8.920180666433653e-16, -8.92478658317932e-17, -8.92478658317932e-17, -8.92478658317932e-17, -1.6929642664672612e-26, -3.881818184146136e-17, -8.92478658317932e-17, 8.922618031428675e-16, 1.669299767318975e-17, 1.669299767318975e-17, -1.6929642664672612e-26, -3.881818184146136e-17, -8.92478658317932e-17, 1.669299767318975e-17, 8.922618031428675e-16, 1.669299767318975e-17, -1.6929642664672612e-26, -3.881818184146136e-17, -8.92478658317932e-17, 1.669299767318975e-17, 1.669299767318975e-17, 8.922618031428675e-16, 
0.07216878418541067, 2.184479142614836e-15, 1.151270483756313e-11, -1.8610120567131596e-11, -1.8610120567131596e-11, -1.8610120567131596e-11, 2.1844791426148358e-15, 4.78506556306281, 0.7413803803373032, 1.8474073036197075e-15, 1.8474073036197075e-15, 1.8474073036197075e-15, 1.1512704837563126e-11, 0.7413803803373029, 1.1189469149298, 5.589230760384249e-16, 5.589230760384249e-16, 5.589230760384249e-16, -1.8610120567131593e-11, 1.8474073036197067e-15, 5.589230760384249e-16, 1.1158138523979704, 1.664587206515389e-30, 1.664587206515389e-30, -1.8610120567131593e-11, 1.8474073036197067e-15, 5.589230760384249e-16, 1.664587206515389e-30, 1.1158138523979704, 1.664587206515389e-30, -1.8610120567131593e-11, 1.8474073036197067e-15, 5.589230760384249e-16, 1.664587206515389e-30, 1.664587206515389e-30, 1.1158138523979704, 
0.01708263539810581, 4.2553798748513794e-16, 2.7250778192552604e-12, -4.405077201779597e-12, -4.405077201779597e-12, -4.405077201779597e-12, 4.25537987485138e-16, 0.7413803803373027, 0.13687339129775566, 3.3490845857242015e-16, 3.3490845857242015e-16, 3.3490845857242015e-16, 2.7250778192552604e-12, 0.1368733912977556, 0.2566334071998988, 1.2730307720343426e-16, 1.2730307720343426e-16, 1.2730307720343426e-16, -4.405077201779597e-12, 3.349084585724201e-16, 1.2730307720343428e-16, 0.2566839963938211, 2.837476178358667e-31, 2.837476178358667e-31, -4.405077201779597e-12, 3.349084585724201e-16, 1.2730307720343428e-16, 2.837476178358667e-31, 0.2566839963938211, 2.837476178358667e-31, -4.405077201779597e-12, 3.349084585724201e-16, 1.2730307720343428e-16, 2.837476178358667e-31, 2.837476178358667e-31, 0.2566839963938211, 
-0.000152732992093364, -3.881818184146136e-17, -7.946365155061487e-14, 1.3582477533571906e-13, 1.150745348691181e-13, 1.150745348691181e-13, -3.881818184146136e-17, 1.8474073036197067e-15, 3.349084585724201e-16, 0.024477413140958703, 1.2016361189427825e-30, 1.2016361189427825e-30, -7.946365155061487e-14, 3.349084585724201e-16, 6.561716020686411e-16, 0.037808608975303734, 6.950646403974642e-31, 6.950646403974642e-31, 1.3582477533571906e-13, 0.024477413140958703, 0.03780860897530374, 6.99140404428425e-16, 3.164311642652069e-17, 3.164311642652069e-17, 1.150745348691181e-13, 1.2016361189427825e-30, 6.950646403974641e-31, 3.164311642652069e-17, 6.358541715753833e-16, 1.184385780501129e-45, 1.150745348691181e-13, 1.2016361189427825e-30, 6.950646403974641e-31, 3.164311642652069e-17, 1.184385780501129e-45, 6.358541715753833e-16, 
-0.000152732992093364, -3.881818184146136e-17, -7.946365155061487e-14, 1.150745348691181e-13, 1.3582477533571906e-13, 1.150745348691181e-13, -3.881818184146136e-17, 1.8474073036197067e-15, 3.349084585724201e-16, 1.2016361189427825e-30, 0.024477413140958703, 1.2016361189427825e-30, -7.946365155061487e-14, 3.349084585724201e-16, 6.561716020686411e-16, 6.950646403974642e-31, 0.037808608975303734, 6.950646403974642e-31, 1.150745348691181e-13, 1.2016361189427825e-30, 6.950646403974641e-31, 6.358541715753833e-16, 3.164311642652069e-17, 1.184385780501129e-45, 1.3582477533571906e-13, 0.024477413140958703, 0.03780860897530374, 3.164311642652069e-17, 6.99140404428425e-16, 3.164311642652069e-17, 1.150745348691181e-13, 1.2016361189427825e-30, 6.950646403974641e-31, 1.184385780501129e-45, 3.164311642652069e-17, 6.358541715753833e-16, 
-0.000152732992093364, -3.881818184146136e-17, -7.946365155061487e-14, 1.150745348691181e-13, 1.150745348691181e-13, 1.3582477533571906e-13, -3.881818184146136e-17, 1.8474073036197067e-15, 3.349084585724201e-16, 1.2016361189427825e-30, 1.2016361189427825e-30, 0.024477413140958703, -7.946365155061487e-14, 3.349084585724201e-16, 6.561716020686411e-16, 6.950646403974642e-31, 6.950646403974642e-31, 0.037808608975303734, 1.150745348691181e-13, 1.2016361189427825e-30, 6.950646403974641e-31, 6.358541715753833e-16, 1.184385780501129e-45, 3.164311642652069e-17, 1.150745348691181e-13, 1.2016361189427825e-30, 6.950646403974641e-31, 1.184385780501129e-45, 6.358541715753833e-16, 3.164311642652069e-17, 1.3582477533571906e-13, 0.024477413140958703, 0.03780860897530374, 3.164311642652069e-17, 3.164311642652069e-17, 6.99140404428425e-16, 
5.1026021939451584e-12, 1.051326389715241e-26, 1.417481549327364e-21, -2.4091087190047952e-21, -2.4091087190047952e-21, -2.4091087190047952e-21, 1.051326389715241e-26, 1.1512704837563123e-11, 2.7250778192552604e-12, -7.946365155061487e-14, -7.946365155061487e-14, -7.946365155061487e-14, 1.417481549327364e-21, 2.7250778192552604e-12, 1.1507004690826273e-11, -9.980982538000196e-13, -9.980982538000196e-13, -9.980982538000196e-13, -2.4091087190047952e-21, -7.946365155061488e-14, -9.980982538000196e-13, 1.1506575748930733e-11, 1.8554495532496332e-13, 1.8554495532496332e-13, -2.4091087190047952e-21, -7.946365155061488e-14, -9.980982538000196e-13, 1.8554495532496332e-13, 1.1506575748930733e-11, 1.8554495532496332e-13, -2.4091087190047952e-21, -7.946365155061488e-14, -9.980982538000196e-13, 1.8554495532496332e-13, 1.8554495532496332e-13, 1.1506575748930733e-11, 
0.017082635398105804, 4.255379874851376e-16, 2.7250778192552596e-12, -4.4050772017795975e-12, -4.4050772017795975e-12, -4.4050772017795975e-12, 4.255379874851377e-16, 0.7413803803373034, 0.1368733912977557, 3.3490845857242015e-16, 3.3490845857242015e-16, 3.3490845857242015e-16, 2.725077819255259e-12, 0.13687339129775564, 0.2566334071998988, 1.2730307720343426e-16, 1.2730307720343426e-16, 1.2730307720343426e-16, -4.4050772017795975e-12, 3.349084585724201e-16, 1.2730307720343428e-16, 0.25668399639382095, 2.837476178358667e-31, 2.837476178358667e-31, -4.4050772017795975e-12, 3.349084585724201e-16, 1.2730307720343428e-16, 2.837476178358667e-31, 0.25668399639382095, 2.837476178358667e-31, -4.4050772017795975e-12, 3.349084585724201e-16, 1.2730307720343428e-16, 2.837476178358667e-31, 2.837476178358667e-31, 0.25668399639382095, 
0.07216878493079162, 8.920180666433649e-16, 1.1507004690826275e-11, -1.8605515587788365e-11, -1.8605515587788365e-11, -1.8605515587788365e-11, 8.920180666433649e-16, 1.1189469149297997, 0.25663340719989874, 6.561716020686412e-16, 6.561716020686412e-16, 6.561716020686412e-16, 1.1507004690826275e-11, 0.2566334071998988, 0.8172063654514495, 0.0, 0.0, 0.0, -1.8605515587788365e-11, 6.561716020686411e-16, 0.0, 0.8170226432280128, 0.0, 0.0, -1.8605515587788365e-11, 6.561716020686411e-16, 0.0, 0.0, 0.8170226432280128, 0.0, -1.8605515587788365e-11, 6.561716020686411e-16, 0.0, 0.0, 0.0, 0.8170226432280128, 
-0.0019280277056418268, -8.924786583179323e-17, -9.9809825380002e-13, 1.708399128358932e-12, 1.4491112343398354e-12, 1.4491112343398354e-12, -8.924786583179324e-17, 5.589230760384248e-16, 1.2730307720343428e-16, 0.037808608975303734, 6.950646403974643e-31, 6.950646403974643e-31, -9.9809825380002e-13, 1.273030772034343e-16, 0.0, 0.18051840048007795, 0.0, 0.0, 1.708399128358932e-12, 0.03780860897530373, 0.18051840048007797, 0.0, 0.0, 0.0, 1.4491112343398356e-12, 6.950646403974643e-31, 0.0, 0.0, 0.0, 0.0, 1.4491112343398356e-12, 6.950646403974643e-31, 0.0, 0.0, 0.0, 0.0, 
-0.0019280277056418268, -8.924786583179323e-17, -9.9809825380002e-13, 1.4491112343398354e-12, 1.708399128358932e-12, 1.4491112343398354e-12, -8.924786583179324e-17, 5.589230760384248e-16, 1.2730307720343428e-16, 6.950646403974643e-31, 0.037808608975303734, 6.950646403974643e-31, -9.9809825380002e-13, 1.273030772034343e-16, 0.0, 0.0, 0.18051840048007795, 0.0, 1.4491112343398356e-12, 6.950646403974643e-31, 0.0, 0.0, 0.0, 0.0, 1.708399128358932e-12, 0.03780860897530373, 0.18051840048007797, 0.0, 0.0, 0.0, 1.4491112343398356e-12, 6.950646403974643e-31, 0.0, 0.0, 0.0, 0.0, 
-0.0019280277056418268, -8.924786583179323e-17, -9.9809825380002e-13, 1.4491112343398354e-12, 1.4491112343398354e-12, 1.708399128358932e-12, -8.924786583179324e-17, 5.589230760384248e-16, 1.2730307720343428e-16, 6.950646403974643e-31, 6.950646403974643e-31, 0.037808608975303734, -9.9809825380002e-13, 1.273030772034343e-16, 0.0, 0.0, 0.0, 0.18051840048007795, 1.4491112343398356e-12, 6.950646403974643e-31, 0.0, 0.0, 0.0, 0.0, 1.4491112343398356e-12, 6.950646403974643e-31, 0.0, 0.0, 0.0, 0.0, 1.708399128358932e-12, 0.03780860897530373, 0.18051840048007797, 0.0, 0.0, 0.0, 
-8.866040348460655e-12, -1.6929642664672612e-26, -2.4091087190047956e-21, 4.196583276467744e-21, 4.0939607446465595e-21, 4.0939607446465595e-21, -1.692964266467261e-26, -1.8610120567131596e-11, -4.405077201779598e-12, 1.3582477533571906e-13, 1.150745348691181e-13, 1.150745348691181e-13, -2.4091087190047956e-21, -4.4050772017795975e-12, -1.8605515587788372e-11, 1.7083991283589316e-12, 1.4491112343398354e-12, 1.4491112343398354e-12, 4.196583276467743e-21, 1.3582477533571906e-13, 1.7083991283589318e-12, -1.866996636040683e-11, -2.891914359470458e-13, -2.891914359470458e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.891914359470458e-13, -1.8572666436284773e-11, -2.405414738860167e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.891914359470458e-13, -2.405414738860167e-13, -1.8572666436284773e-11, 
-0.00015273299209336394, -3.8818181841461365e-17, -7.946365155061485e-14, 1.3582477533571908e-13, 1.150745348691181e-13, 1.150745348691181e-13, -3.881818184146136e-17, 1.8474073036197067e-15, 3.3490845857242015e-16, 0.024477413140958706, 1.2016361189427827e-30, 1.2016361189427827e-30, -7.946365155061485e-14, 3.3490845857242015e-16, 6.561716020686412e-16, 0.03780860897530376, 6.950646403974642e-31, 6.950646403974642e-31, 1.3582477533571908e-13, 0.024477413140958706, 0.03780860897530377, 6.99140404428425e-16, 3.1643116426520694e-17, 3.1643116426520694e-17, 1.150745348691181e-13, 1.2016361189427827e-30, 6.950646403974641e-31, 3.1643116426520694e-17, 6.358541715753835e-16, 1.1843857805011287e-45, 1.150745348691181e-13, 1.2016361189427827e-30, 6.950646403974641e-31, 3.1643116426520694e-17, 1.1843857805011287e-45, 6.358541715753835e-16, 
-0.0019280277056418268, -8.92478658317932e-17, -9.9809825380002e-13, 1.708399128358932e-12, 1.4491112343398354e-12, 1.4491112343398354e-12, -8.924786583179322e-17, 5.58923076038425e-16, 1.2730307720343428e-16, 0.03780860897530373, 6.950646403974641e-31, 6.950646403974641e-31, -9.980982538000198e-13, 1.273030772034343e-16, 0.0, 0.1805184004800779, 0.0, 0.0, 1.7083991283589322e-12, 0.03780860897530372, 0.18051840048007786, 0.0, 0.0, 0.0, 1.4491112343398354e-12, 6.950646403974642e-31, 0.0, 0.0, 0.0, 0.0, 1.4491112343398354e-12, 6.950646403974642e-31, 0.0, 0.0, 0.0, 0.0, 
0.07216878440005933, 8.922618031428673e-16, 1.1506575748930733e-11, -1.8669966360406834e-11, -1.857266643628477e-11, -1.857266643628477e-11, 8.922618031428671e-16, 1.1158138523979704, 0.2566839963938211, 6.991404044284248e-16, 6.358541715753834e-16, 6.358541715753834e-16, 1.1506575748930732e-11, 0.25668399639382106, 0.817022643228013, 0.0, 0.0, 0.0, -1.8669966360406837e-11, 6.991404044284249e-16, 0.0, 0.8801591277387371, 0.0, 0.0, -1.857266643628477e-11, 6.358541715753834e-16, 0.0, 0.0, 0.7852702337972611, 0.0, -1.857266643628477e-11, 6.358541715753834e-16, 0.0, 0.0, 0.0, 0.7852702337972611, 
0.00011167016843450897, 1.669299767318975e-17, 1.8554495532496332e-13, -2.8919143594704597e-13, -2.8919143594704597e-13, -2.4054147388601677e-13, 1.6692997673189746e-17, 1.6645872065153893e-30, 2.8374761783586664e-31, 3.164311642652069e-17, 3.164311642652069e-17, 1.184385780501129e-45, 1.8554495532496332e-13, 2.837476178358667e-31, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.04744444697073826, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.04744444697073826, 0.0, 0.0, -2.4054147388601677e-13, 1.1843857805011292e-45, 0.0, 0.0, 0.0, 0.0, 
0.00011167016843450897, 1.669299767318975e-17, 1.8554495532496332e-13, -2.8919143594704597e-13, -2.405414738860168e-13, -2.8919143594704597e-13, 1.6692997673189746e-17, 1.6645872065153893e-30, 2.8374761783586664e-31, 3.164311642652069e-17, 1.184385780501129e-45, 3.164311642652069e-17, 1.8554495532496332e-13, 2.837476178358667e-31, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.0, 0.04744444697073826, -2.405414738860168e-13, 1.1843857805011292e-45, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.04744444697073826, 0.0, 0.0, 
-8.866040348460655e-12, -1.6929642664672612e-26, -2.4091087190047956e-21, 4.0939607446465595e-21, 4.196583276467744e-21, 4.0939607446465595e-21, -1.692964266467261e-26, -1.8610120567131596e-11, -4.405077201779598e-12, 1.150745348691181e-13, 1.3582477533571906e-13, 1.150745348691181e-13, -2.4091087190047956e-21, -4.4050772017795975e-12, -1.8605515587788372e-11, 1.4491112343398354e-12, 1.7083991283589316e-12, 1.4491112343398354e-12, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -1.857266643628477e-11, -2.891914359470458e-13, -2.405414738860167e-13, 4.196583276467743e-21, 1.3582477533571906e-13, 1.7083991283589318e-12, -2.891914359470458e-13, -1.866996636040683e-11, -2.891914359470458e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.405414738860167e-13, -2.891914359470458e-13, -1.8572666436284773e-11, 
-0.00015273299209336394, -3.8818181841461365e-17, -7.946365155061485e-14, 1.150745348691181e-13, 1.3582477533571908e-13, 1.150745348691181e-13, -3.881818184146136e-17, 1.8474073036197067e-15, 3.3490845857242015e-16, 1.2016361189427827e-30, 0.024477413140958706, 1.2016361189427827e-30, -7.946365155061485e-14, 3.3490845857242015e-16, 6.561716020686412e-16, 6.950646403974642e-31, 0.03780860897530376, 6.950646403974642e-31, 1.150745348691181e-13, 1.2016361189427827e-30, 6.950646403974641e-31, 6.358541715753835e-16, 3.1643116426520694e-17, 1.1843857805011287e-45, 1.3582477533571908e-13, 0.024477413140958706, 0.03780860897530377, 3.1643116426520694e-17, 6.99140404428425e-16, 3.1643116426520694e-17, 1.150745348691181e-13, 1.2016361189427827e-30, 6.950646403974641e-31, 1.1843857805011287e-45, 3.1643116426520694e-17, 6.358541715753835e-16, 
-0.0019280277056418268, -8.92478658317932e-17, -9.9809825380002e-13, 1.4491112343398354e-12, 1.708399128358932e-12, 1.4491112343398354e-12, -8.924786583179322e-17, 5.58923076038425e-16, 1.2730307720343428e-16, 6.950646403974641e-31, 0.03780860897530373, 6.950646403974641e-31, -9.980982538000198e-13, 1.273030772034343e-16, 0.0, 0.0, 0.1805184004800779, 0.0, 1.4491112343398354e-12, 6.950646403974642e-31, 0.0, 0.0, 0.0, 0.0, 1.7083991283589322e-12, 0.03780860897530372, 0.18051840048007786, 0.0, 0.0, 0.0, 1.4491112343398354e-12, 6.950646403974642e-31, 0.0, 0.0, 0.0, 0.0, 
0.00011167016843450897, 1.669299767318975e-17, 1.8554495532496332e-13, -2.8919143594704597e-13, -2.8919143594704597e-13, -2.4054147388601677e-13, 1.6692997673189746e-17, 1.6645872065153893e-30, 2.8374761783586664e-31, 3.164311642652069e-17, 3.164311642652069e-17, 1.184385780501129e-45, 1.8554495532496332e-13, 2.837476178358667e-31, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.04744444697073826, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.04744444697073826, 0.0, 0.0, -2.4054147388601677e-13, 1.1843857805011292e-45, 0.0, 0.0, 0.0, 0.0, 
0.07216878440005933, 8.922618031428673e-16, 1.1506575748930733e-11, -1.857266643628477e-11, -1.8669966360406834e-11, -1.857266643628477e-11, 8.922618031428671e-16, 1.1158138523979704, 0.2566839963938211, 6.358541715753834e-16, 6.991404044284248e-16, 6.358541715753834e-16, 1.1506575748930732e-11, 0.25668399639382106, 0.817022643228013, 0.0, 0.0, 0.0, -1.8572666436284773e-11, 6.358541715753834e-16, 0.0, 0.7852702337972611, 0.0, 0.0, -1.8669966360406837e-11, 6.991404044284249e-16, 0.0, 0.0, 0.8801591277387371, 0.0, -1.857266643628477e-11, 6.358541715753834e-16, 0.0, 0.0, 0.0, 0.7852702337972611, 
0.00011167016843450897, 1.669299767318975e-17, 1.8554495532496332e-13, -2.405414738860168e-13, -2.8919143594704597e-13, -2.8919143594704597e-13, 1.6692997673189746e-17, 1.6645872065153893e-30, 2.8374761783586664e-31, 1.184385780501129e-45, 3.164311642652069e-17, 3.164311642652069e-17, 1.8554495532496332e-13, 2.837476178358667e-31, 0.0, 0.0, 0.0, 0.0, -2.405414738860168e-13, 1.1843857805011292e-45, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.0, 0.04744444697073826, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.04744444697073826, 0.0, 
-8.866040348460655e-12, -1.6929642664672612e-26, -2.4091087190047956e-21, 4.0939607446465595e-21, 4.0939607446465595e-21, 4.196583276467744e-21, -1.692964266467261e-26, -1.8610120567131596e-11, -4.405077201779598e-12, 1.150745348691181e-13, 1.150745348691181e-13, 1.3582477533571906e-13, -2.4091087190047956e-21, -4.4050772017795975e-12, -1.8605515587788372e-11, 1.4491112343398354e-12, 1.4491112343398354e-12, 1.7083991283589316e-12, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -1.857266643628477e-11, -2.405414738860167e-13, -2.891914359470458e-13, 4.093960744646559e-21, 1.150745348691181e-13, 1.4491112343398352e-12, -2.405414738860167e-13, -1.857266643628477e-11, -2.891914359470458e-13, 4.196583276467743e-21, 1.3582477533571906e-13, 1.7083991283589318e-12, -2.891914359470458e-13, -2.891914359470458e-13, -1.866996636040683e-11, 
-0.00015273299209336394, -3.8818181841461365e-17, -7.946365155061485e-14, 1.150745348691181e-13, 1.150745348691181e-13, 1.3582477533571908e-13, -3.881818184146136e-17, 1.8474073036197067e-15, 3.3490845857242015e-16, 1.2016361189427827e-30, 1.2016361189427827e-30, 0.024477413140958706, -7.946365155061485e-14, 3.3490845857242015e-16, 6.561716020686412e-16, 6.950646403974642e-31, 6.950646403974642e-31, 0.03780860897530376, 1.150745348691181e-13, 1.2016361189427827e-30, 6.950646403974641e-31, 6.358541715753835e-16, 1.1843857805011287e-45, 3.1643116426520694e-17, 1.150745348691181e-13, 1.2016361189427827e-30, 6.950646403974641e-31, 1.1843857805011287e-45, 6.358541715753835e-16, 3.1643116426520694e-17, 1.3582477533571908e-13, 0.024477413140958706, 0.03780860897530377, 3.1643116426520694e-17, 3.1643116426520694e-17, 6.99140404428425e-16, 
-0.0019280277056418268, -8.92478658317932e-17, -9.9809825380002e-13, 1.4491112343398354e-12, 1.4491112343398354e-12, 1.708399128358932e-12, -8.924786583179322e-17, 5.58923076038425e-16, 1.2730307720343428e-16, 6.950646403974641e-31, 6.950646403974641e-31, 0.03780860897530373, -9.980982538000198e-13, 1.273030772034343e-16, 0.0, 0.0, 0.0, 0.1805184004800779, 1.4491112343398354e-12, 6.950646403974642e-31, 0.0, 0.0, 0.0, 0.0, 1.4491112343398354e-12, 6.950646403974642e-31, 0.0, 0.0, 0.0, 0.0, 1.7083991283589322e-12, 0.03780860897530372, 0.18051840048007786, 0.0, 0.0, 0.0, 
0.00011167016843450897, 1.669299767318975e-17, 1.8554495532496332e-13, -2.8919143594704597e-13, -2.405414738860168e-13, -2.8919143594704597e-13, 1.6692997673189746e-17, 1.6645872065153893e-30, 2.8374761783586664e-31, 3.164311642652069e-17, 1.184385780501129e-45, 3.164311642652069e-17, 1.8554495532496332e-13, 2.837476178358667e-31, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.0, 0.04744444697073826, -2.405414738860168e-13, 1.1843857805011292e-45, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.04744444697073826, 0.0, 0.0, 
0.00011167016843450897, 1.669299767318975e-17, 1.8554495532496332e-13, -2.405414738860168e-13, -2.8919143594704597e-13, -2.8919143594704597e-13, 1.6692997673189746e-17, 1.6645872065153893e-30, 2.8374761783586664e-31, 1.184385780501129e-45, 3.164311642652069e-17, 3.164311642652069e-17, 1.8554495532496332e-13, 2.837476178358667e-31, 0.0, 0.0, 0.0, 0.0, -2.405414738860168e-13, 1.1843857805011292e-45, 0.0, 0.0, 0.0, 0.0, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.0, 0.04744444697073826, -2.8919143594704597e-13, 3.1643116426520694e-17, 0.0, 0.0, 0.04744444697073826, 0.0, 
0.07216878440005933, 8.922618031428673e-16, 1.1506575748930733e-11, -1.857266643628477e-11, -1.857266643628477e-11, -1.8669966360406834e-11, 8.922618031428671e-16, 1.1158138523979704, 0.2566839963938211, 6.358541715753834e-16, 6.358541715753834e-16, 6.991404044284248e-16, 1.1506575748930732e-11, 0.25668399639382106, 0.817022643228013, 0.0, 0.0, 0.0, -1.8572666436284773e-11, 6.358541715753834e-16, 0.0, 0.7852702337972611, 0.0, 0.0, -1.8572666436284773e-11, 6.358541715753834e-16, 0.0, 0.0, 0.7852702337972611, 0.0, -1.8669966360406837e-11, 6.991404044284249e-16, 0.0, 0.0, 0.0, 0.8801591277387371};


float LU_test[cgto_number * cgto_number] = {
    1.f, 0.f, 0.f, 0.f, 0.f,0.f,
    0.f, 1.f, 0.236704f, 0.f, 0.f,0.f,
    0.f, 0.236704f  , 0.94397122f, 0.f, 0.f,0.f,
    0.f, 0.f, 0.f, 1.f, 0.f,0.f,
    0.f, 0.f, 0.f, 0.f, 1.f,0.f,
    0.f, 0.f, 0.f, 0.f, 0.f,1.f};
float INV_test[ cgto_number*cgto_number] = { 
    1.f, -0.f, -0.f, -0.f, 0.f,-0.f,
    0.f,  1.05935433f, -0.25075341f, -0.f , 0.f,-0.f,
    0.f, -0.25075341f,  1.05935433f,  0.f , 0.f,-0.f,
    0.f, 0.f,  0.f,  1.f , 0.f,-0.f,
    0.f, 0.f,  0.f,  0.f , 1.f,-0.f,
    0.f, 0.f,  0.f,  0.f , 0.f,1.0f};

float Eval_test[cgto_number] = {0.763297,1.000000,1.000000,1.000000,1.000000,1.236704};
float Evec_test[cgto_number * cgto_number] = {
    0.000000,1.000000,0.000000,0.000000,0.000000,0.000000,
    -0.707107,0.000000,0.000000,0.000000,0.000000,0.707107,
    0.707107,0.000000,0.000000,0.000000,0.000000,0.707107,
    0.000000,0.000000,1.000000,0.000000,0.000000,0.000000,
    0.000000,0.000000,0.000000,1.000000,0.000000,0.000000,
    0.000000,0.000000,0.000000,0.000000,1.000000,0.000000
    };