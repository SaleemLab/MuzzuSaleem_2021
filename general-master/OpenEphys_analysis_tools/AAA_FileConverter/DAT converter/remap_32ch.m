%% Tomaso Muzzu - UCL - 28 Feb 2018

% function to reorder the channels following the order:
% the saved 32 channels from the OpenEphys software are re-mapped to the
% physical channel positions on the shanks of the probe used. Each channel
% saved by the OpenEphys software as CH*.continuous will be mapped through 
% the OM32 connector pins of the RHD2132 amplifier board, the OM32-A32 
% connector, the A32 connector of the probe and, finally, the map of the
% channels on the shanks of the probe.
% ProbeBrand=1(NeuroNexus), ProbeBrand=2(CambridgeNeurotech)

function Ch_map = remap_32ch(ProbeBrand)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANNELS' ORDER ON A32 CONNECTOR
switch ProbeBrand
    case 1
        Connector =  [32  0	0	11 ...
            30	0	0	9 ...
            31	0	0	7 ...
            28	0	0	5 ...
            29	26	1	3 ...
            27	24	4	2 ...
            25	20	13	6 ...
            22	19	14	8 ...
            23	18	15	10 ...
            21	17	16	12];
    case 2
        Connector =  [32  0	0	11 ...
            30	0	0	9 ...
            31	0	0	7 ...
            28	0	0	5 ...
            29	26	1	3 ...
            27	24	4	2 ...
            25	20	13	6 ...
            22	19	14	8 ...
            23	18	15	10 ...
            21	17	16	12];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A32-OM32 ADAPTOR MAP
% A32 female connector (ordered as they connect to the male connector above)
A32_OM32_adapter = ...
    [32	0	0	1 ...
    31	0	0	2 ...
    30	0	0	3 ...
    29	0	0	4 ...
    28	17	16	5 ...
    27	18	15	6 ...
    26	19	14	7 ...
    25	20	13	8 ...
    24	21	12	9 ...
    23	22	11	10];

% OM32 connector
% order from left to right with 'NEURONEXUS' sign on top (see figure for ref.)
OM32_connector = ...
[	23	25	27	29	31	19	17	21	11	15	13	1	3	5	7	9	 ...
	24	26	28	30	32	20	18	22	12	16	14	2	4	6	8	10	];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RHD2132 amplifier board
% ordered as they connect to the female OM32 connector above
RHD2132 = ...
[	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	 ...
	8	7	6	5	4	3	2	1	32	31	30	29	28	27	26	25	];

%% reorder channels now
% OpenEphys channel 1 signal comes from the 24th pin of the RHD2132
% connector. This pin is connected to the 24th adapter pin that is carrying 
% the 22nd channel signal. The signal of channel 22 is going through the
% 38th pin of the A32 connector of the probe that carries the 17th channel
% signal. The 17th channel in on the third shank. Let's translate this into
% matlab language, considering channel 1:
% find(A32_connector==1) 
% A32_OM32_adapter(find(A32_connector==1))
% find(OM32_connector==A32_OM32_adapter(find(A32_connector==1)))
% RHD2132(find(OM32_connector==A32_OM32_adapter(find(A32_connector==1))))
for i = 1:32
    Ch_map(i) = RHD2132(find(OM32_connector==A32_OM32_adapter(find(Connector==i))));
end
% [~,Ch_map] = sort(Ch_map);

% % try the order
% A1 = Ephys_data(1:100,:);
% A2 = A1(:,Ch_map);
% alt1Aman = [18    26    20    30    28    27    22    21    23    25 ...
% 29    31    32    19    24    17    15     1    13     9     5     2 ...
% 22    21    23    25    29    31    32    19    24    20];

% alt2Aman = [15     1    13     9     5     2     8    12    10     6 ...
% 3     7    11    14     4    16    18    26    31    30    28    27 ...
% 22    21    23    25    29    31    32    19    24    20]; % This should be
% correct

end
