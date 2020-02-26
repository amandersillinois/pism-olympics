CDF       
      time  �         Comments      }Revised on 02 February 2010 by Glen Granzow at the University of Montana to include climate data from J. Ettema et al. (2009)      Conventions       CF-1.3     Creators      <Jesse Johnson, Brian Hand, Tim Bocek - University of Montana   History      �Wed Feb 26 14:20:07 2020: ncatted -O -a units,delta_T,m,c,Kelvin pism_dT.nc
Wed Feb 26 14:20:07 2020: ncatted -O -a calendar,time,c,c,365_day pism_dT.nc
Wed Feb 26 14:20:07 2020: ncatted -O -a units,time,m,c,years since 1-1-1 pism_dT.nc
Wed Feb 26 14:20:06 2020: ncap2 -O -s time=-time pism_dT.nc pism_dT.nc
Wed Feb 26 14:20:06 2020: ncpdq -O --rdr=-time pism_dT.nc pism_dT.nc
Wed Feb 26 14:20:06 2020: ncrename -O -v oisotopestimes,time pism_dT.nc
Wed Feb 26 14:20:06 2020: ncrename -O -v temp_time_series,delta_T pism_dT.nc
Wed Feb 26 14:20:06 2020: ncrename -O -d oisotopestimes,time pism_dT.nc
Wed Feb 26 14:20:06 2020: ncks -O -v oisotopestimes,temp_time_series Greenland_5km_v1.1.nc pism_dT.nc
Original data set created February 2009    Title         Greenland Standard Data Set    NCO       `netCDF Operators version 4.9.2 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)          time                	long_name         Time (years before present)    standard_name         time   units         years since 1-1-1      calendar      365_day      �  l   delta_T                 	reference        �[1] Johnsen, S.J., H.B. Clausen, W. Dansgaard, N.S. Gundestrup, C.U. Hammer, U. Andersen, K.K. Andersen, C.S. Hvidberg, D. Dahl-Jensen, J.P. Steffensen, H. Shoji, A.E. Sveinbj-rnsdUttir, J.W.C. White, J. Jouzel, and D. Fisher (1997), "The d18O Record Along the Greenland Ice Core Project Deep Ice Core and the Problem of Possible Eemian Climatic Instability"; Journal of Geophysical Research 102:26397-26410.  [2] Dansgaard, W., S.J. Johnsen, H.B. Clausen, D. Dahl-Jensen, N.S. Gundestrup, C.U. Hammer, C.S. Hvidberg, J.P. Steffensen, A.E. Sveinbj-rnsdUttir, J. Jouzel, and G.C. Bond (1993), "Evidence for General Instability of Past Climate from a 250 kyr Ice-core Record"; Nature 264:218-220.  [3] GRIP Members (1993), "Climate Instability During the Last Interglacial Period Recorded in the GRIP Ice Core"; Nature 364:203-207.  [4] Grootes, P.M., M. Stuiver, J.W.C. White, S.J. Johnsen, and J. Jouzel (1993), "Comparison of Oxygen Isotope Records from the GISP2 and GRIP Greenland Ice Cores"; Nature 366:552-554.  [5] Dansgaard, W., J.W.C. White, and S.J. Johnsen (1989), "The Abrupt Termination of the Younger Dryas Climate Event"; Nature 339:532-533.    computed_by_formula       i2.4(oisotope_time_series+34.83) (except for t=0 at which the temperature variation is zero by definition)      	long_name         $Temperature (variation from present)   standard_name         !land_ice_temperature_at_firn_base      units         Kelvin       �  ���$ ��� ��� �� ��\ ��* ��� ��� �� ��b ��0 ��� ��� �� ��h ��6 �� ��� �� ��n ��< ��
 ��� �� ��t ��B �� ��� �� ��z ��H �� ��� ��� �� ��N �� ��� �� �� ��T ��" ��� �� �� ��Z ��( ��� ��� �� ��` ��. ��� ��� �� ��f ��4 �� ��� �� ��l ��: �� ��� �� ��r ��@ �� ��� �� ��x ��F �� ��� �� ��~ ��L �� ��� �� �� ��R ��  ��� �� �� ��X ��& ��� ��� �� ��^ ��, ��� ��� �� ��d ��2 ��  ��� ��� ��j ��8 �� ��� �ߢ ��p ��> �� ��� �ި ��v ��D �� ��� �ݮ ��| ��J �� ��� �ܴ �܂ ��P �� ��� �ۺ �ۈ ��V ��$ ��� ��� �ڎ ��\ ��* ��� ��� �ٔ ��b ��0 ��� ��� �ؚ ��h ��6 �� ��� �נ ��n ��< ��
 ��� �֦ ��t ��B �� ��� �լ ��z ��H �� ��� �Բ �Ԁ ��N �� ��� �Ӹ �ӆ ��T ��" ��� �Ҿ �Ҍ ��Z ��( ��� ��� �ђ ��` ��. ��� ��� �И ��f ��4 �� ��� �Ϟ ��l ��: �� ��� �Τ ��r ��@ �� ��� �ͪ ��x ��F �� ��� �̰ ��~ ��L �� ��� �˶ �˄ ��R ��  ��� �ʼ �ʊ ��X ��& ��� ��� �ɐ ��^ ��, ��� ��� �Ȗ ��d ��2 ��  ��� �ǜ ��j ��8 �� ��� �Ƣ ��p ��> �� ��� �Ũ ��v ��D �� ��� �Į ��| ��J �� ��� �ô �Â ��P �� ��� �º � ��V ��$ ��� ��� ��� ��\ ��* ��� ��� ��� ��b ��0 ǿ� ǿ� ǿ� ǿh ǿ6 ǿ Ǿ� Ǿ� Ǿn Ǿ< Ǿ
 ǽ� ǽ� ǽt ǽB ǽ Ǽ� Ǽ� Ǽz ǼH Ǽ ǻ� ǻ� ǻ� ǻN ǻ Ǻ� Ǻ� Ǻ� ǺT Ǻ" ǹ� ǹ� ǹ� ǹZ ǹ( Ǹ� Ǹ� Ǹ� Ǹ` Ǹ. Ƿ� Ƿ� Ƿ� Ƿf Ƿ4 Ƿ Ƕ� Ƕ� Ƕl Ƕ: Ƕ ǵ� ǵ� ǵr ǵ@ ǵ Ǵ� Ǵ� Ǵx ǴF Ǵ ǳ� ǳ� ǳ~ ǳL ǳ ǲ� ǲ� ǲ� ǲR ǲ  Ǳ� Ǳ� Ǳ� ǱX Ǳ& ǰ� ǰ� ǰ� ǰ^ ǰ, ǯ� ǯ� ǯ� ǯd ǯ2 ǯ  Ǯ� Ǯ� Ǯj Ǯ8 Ǯ ǭ� ǭ� ǭp ǭ> ǭ Ǭ� Ǭ� Ǭv ǬD Ǭ ǫ� ǫ� ǫ| ǫJ ǫ Ǫ� Ǫ� Ǫ� ǪP Ǫ ǩ� ǩ� ǩ� ǩV ǩ$ Ǩ� Ǩ� Ǩ� Ǩ\ Ǩ* ǧ� ǧ� ǧ� ǧb ǧ0 Ǧ� Ǧ� Ǧ� Ǧh Ǧ6 Ǧ ǥ� ǥ� ǥn ǥ< ǥ
 Ǥ� Ǥ� Ǥt ǤB Ǥ ǣ� ǣ� ǣz ǣH ǣ Ǣ� Ǣ� Ǣ� ǢN Ǣ ǡ� ǡ� ǡ� ǡT ǡ" Ǡ� Ǡ� Ǡ� ǠZ Ǡ( ǟ� ǟ� ǟ� ǟ` ǟ. Ǟ� Ǟ� Ǟ� Ǟf Ǟ4 Ǟ ǝ� ǝ� ǝl ǝ: ǝ ǜ� ǜ� ǜr ǜ@ ǜ Ǜ� Ǜ� Ǜx ǛF Ǜ ǚ� ǚ� ǚ~ ǚL ǚ Ǚ� Ǚ� Ǚ� ǙR Ǚ  ǘ� ǘ� ǘ� ǘX ǘ& Ǘ� Ǘ� Ǘ� Ǘ^ Ǘ, ǖ� ǖ� ǖ� ǖd ǖ2 ǖ  Ǖ� Ǖ� Ǖj Ǖ8 Ǖ ǔ� ǔ� ǔp ǔ> ǔ Ǔ� Ǔ� Ǔv ǓD Ǔ ǒ� ǒ� ǒ| ǒJ ǒ Ǒ� Ǒ� Ǒ� ǑP Ǒ ǐ� ǐ� ǐ� ǐV ǐ$ Ǐ� Ǐ� Ǐ� Ǐ\ Ǐ* ǎ� ǎ� ǎ� ǎb ǎ0 Ǎ� Ǎ� Ǎ� Ǎh Ǎ6 Ǎ ǌ� ǌ� ǌn ǌ< ǌ
 ǋ� ǋ� ǋt ǋB ǋ Ǌ� Ǌ� Ǌz ǊH Ǌ ǉ� ǉ� ǉ� ǉN ǉ ǈ� ǈ� ǈ� ǈT ǈ" Ǉ� Ǉ� Ǉ� ǇZ Ǉ( ǆ� ǆ� ǆ� ǆ` ǆ. ǅ� ǅ� ǅ� ǅf ǅ4 ǅ Ǆ� Ǆ� Ǆl Ǆ: Ǆ ǃ� ǃ� ǃr ǃ@ ǃ ǂ� ǂ� ǂx ǂF ǂ ǁ� ǁ� ǁ~ ǁL ǁ ǀ� ǀ� ǀ� ǀR ǀ  �� �x � �~� �~L �}� �}� �}  �|� �|X �{� �{� �{, �z� �zd �z  �y� �y8 �x� �xp �x �w� �wD �v� �v| �v �u� �uP �t� �t� �t$ �s� �s\ �r� �r� �r0 �q� �qh �q �p� �p< �o� �ot �o �n� �nH �m� �m� �m �l� �lT �k� �k� �k( �j� �j` �i� �i� �i4 �h� �hl �h �g� �g@ �f� �fx �f �e� �eL �d� �d� �d  �c� �cX �b� �b� �b, �a� �ad �a  �`� �`8 �_� �_p �_ �^� �^D �]� �]| �] �\� �\P �[� �[� �[$ �Z� �Z\ �Y� �Y� �Y0 �X� �Xh �X �W� �W< �V� �Vt �V �U� �UH �T� �T� �T �S� �ST �R� �R� �R( �Q� �Q` �P� �P� �P4 �O� �Ol �O �N� �N@ �M� �Mx �M �L� �LL �K� �K� �K  �J� �JX �I� �I� �I, �H� �Hd �H  �G� �G8 �F� �Fp �F �E� �ED �D� �D| �D �C� �CP �B� �B� �B$ �A� �A\ �@� �@� �@0 �?� �?h �? �>� �>< �=� �=t �= �<� �<H �;� �;� �; �:� �:T �9� �9� �9( �8� �8` �7� �7� �74 �6� �6l �6 �5� �5@ �4� �4x �4 �3� �3L �2� �2� �2  �1� �1X �0� �0� �0, �/� �/d �/  �.� �.8 �-� �-p �- �,� �,D �+� �+| �+ �*� �*P �)� �)� �)$ �(� �(\ �'� �'� �'0 �&� �&h �& �%� �%< �$� �$t �$ �#� �#H �"� �"� �" �!� �!T � � � � � ( �� �` �� �� �4 �� �l � �� �@ �� �x � �� �L �� �� �  �� �X �� �� �, �� �d �  �� �8 �� �p � �� �D �� �| � �� �P �� �� �$ �� �\ �� �� �0 �� �h � �� �< �� �t � �
� �
H �	� �	� �	 �� �T �� �� �( �� �` �� �� �4 �� �l � �� �@ �� �x � �� �L � � � � �   ��x ��� ��� ��  ��X ��� ��� ��  ��8 ��p ��� ��� �� ��P �� ��� ��� ��0 ��h �� ��� �� ��H �� �� ��� ��( ��` �� ��� �� ��@ ��x �� ��� ��  ��X �� ��� ��  ��8 ��p �ި ��� �� ��P �ۈ ��� ��� ��0 ��h �נ ��� �� ��H �Ԁ �Ӹ ��� ��( ��` �И ��� �� ��@ ��x �̰ ��� ��  ��X �ɐ ��� ��  ��8 ��p �Ũ ��� �� ��P � ��� ��� ��0 ƿh ƾ� ƽ� ƽ ƼH ƻ� ƺ� ƹ� ƹ( Ƹ` Ʒ� ƶ� ƶ Ƶ@ ƴx Ƴ� Ʋ� Ʋ  ƱX ư� Ư� Ư  Ʈ8 ƭp Ƭ� ƫ� ƫ ƪP Ʃ� ƨ� Ƨ� Ƨ0 Ʀh ƥ� Ƥ� Ƥ ƣH Ƣ� ơ� Ơ� Ơ( Ɵ` ƞ� Ɲ� Ɲ Ɯ@ ƛx ƚ� ƙ� ƙ  ƘX Ɨ� Ɩ� Ɩ  ƕ8 Ɣp Ɠ� ƒ� ƒ ƑP Ɛ� Ə� Ǝ� Ǝ0 ƍh ƌ� Ƌ� Ƌ ƊH Ɖ� ƈ� Ƈ� Ƈ( Ɔ` ƅ� Ƅ� Ƅ ƃ@ Ƃx Ɓ� ƀ� ƀ  �~� �}  �{� �z  �xp �v� �uP �s� �r0 �p� �o �m� �k� �j` �h� �g@ �e� �d  �b� �a  �_p �]� �\P �Z� �Y0 �W� �V �T� �R� �Q` �O� �N@ �L� �K  �I� �H  �Fp �D� �CP �A� �@0 �>� �= �;� �9� �8` �6� �5@ �3� �2  �0� �/  �-p �+� �*P �(� �'0 �%� �$ �"� � � �` �� �@ �� �  �� �  �p �� �P �� �0 �� � �	� �� �` �� �@ �� �   ��  ��  ��� ��� �� �� ��` ��@ ��  ��  ��� ��� �נ �Ԁ ��` ��@ ��  ��  ��� ��� ž� Ż� Ÿ` ŵ@ Ų  ů  ū� Ũ� ť� Ţ� ş` Ŝ@ ř  Ŗ  Œ� ŏ� Ō� ŉ� ņ` Ń@ ŀ  �z  �s� �m� �g@ �a  �Z� �T� �N@ �H  �A� �;� �5@ �/  �(� �"� �@ �  �� �	� �@ ��  �� ��  �Ԁ ��  Ļ� į  Ģ� Ė  ĉ� �z  �a  �H  �/  �  ��  ��  Ö  �H  ��  �   ���(��ͨ���,�������&���f�dӧ�P ��i�^��E!���~��J��O���fz��V����5��k%������n���d����S���k�����l$�r�^�q���^�,�<m�>=.�V�`������.���eW������~��A�p��]�I�bw��/-߿�΀���`�]�d����B�?cS�?ZW\=�E�aގ�j�¿ `�?\�?@��@�"�@��@���@���@��~@֤�@�:@Q_@�LK@͕f@� ۿ���(���E��7�ވP�������7�����������\������'�t�����ƿu,��=�j�R���?�N���6b��Ĝ�uS��n4���>y��?��6@�_@�\n@�uW<;>�@�{@Ǯ@��?��K?Ƚ�?��?����>�4?�hs@yD@L3x@PRa?�?�Q�?��=W
=�f�����#���!������ ��ډS��m3��Z���L��5?����(�C����Yi�;��%���!Nv����H-��������z������������ܪ�%���.���&���\)��0�������[���������R�%3��+^e�(�� ��0���F+�Q��G��������.���~~�Ϟ��钘�����ƨ��qS��_�/9�����O�߫���Rt��5���pX����瑖���!����bl��qv��ff��]������;��I����������Z6��+@��H��#��V���v�Z �$�������
�o������{P��\)�~�#Ґ�/��e��.���_����p�j� �b�(��q�����<�+�7���z����<^��P��BW���At�4q&�FS �L&�6����}��oH������*1��~���C ���;��-n���Z������[����������9�����Zr���_��g�����=��������ʊ���+���r��� ����������{��ۥ�������?}��3���j�����+�������6���"z��&K��b����	�#�_���/��/�
	��~�Q�	ӣ�	�B�Z%������������#'���b�u%�.��0R�����m!����Ü������Sq��w��������6��G�1'��8�=q�'H�����-���)�l�%�e�(	7�.V�)�^�*ܙ�*J��'8�'\L�$/��(���45��;M� �S�0L��:�<%�4��1���/y��3W�;}��F��3���3�5V�.+��(���*���7���92��,�}�.9`�7K��7�s�5���8��G���A0!�5�}�1�n�4J� ��$�$�#ޔ�+;S�39��0���$|N�)U\��;�"��ڐ�&$K�m#��Ƴ����٬����̭#������N����������
�������	��f��P�����H�z��R��-�$x�)@D�=�ς�N�R�Z3��a4��gq:�\x��R@�J$��K�_�^ ��B�o�B�!�I��C[�0G��7,��7��(���3Y�2\��"�!Tv�A�#��Ɇ�"���],��K��~;�%��!k�5V��l��8��8��������"��}k�YH��P�������A��U��Y<���������������:���)�[��������G����������F��Ŗ���������A�������m��C����"�������h���d���\�¼��ȮD�Ŵ��ي�������a�◍�����k���Y��\��8�f�����:b�TC��+�]C�W�'ڍ�$cy�)E��8\a�%]�ｔ�����e��+���Z;���� �(�%�,� V�>2f�9Җ�?xi����C�"Pn�?>A�^2��x�	���?��j\�������M��<��V7���_��f���=���������
��'e��7h��ɩ��j���S���-��5v����R�������{���m�����-��Y�̌� ��['�ͥ��-�"e��'7h�)|\�$�j�Z�^��R����T��O:��e������� ���H�������������ݷ��>����~���^���9��3��������E�
��7���o����5q��1'��w=��?��ꄎ��S��ހ%� ��]�H��	���!P��<��39v�W�o�SV�p����e}��֡������F������@��;�{����G5��pY/������5Z�f1�u4��p�(�{�l��'�������՚������P��~<����+�{e������=u��<���-e�t�m�����,��n��w
�����f�w�{����!��{�a�|��~%3�~��������f��u'���������̰��F��Ӭ�p���\���q�7�z��������I��c����"��/��}�n�B���{���\��P{�\1��b��[�?�Y5$�gl��k̔���~��1'��R}�{?������/����r� ��Q�^V�bf��x�Z���1��K,�|Y^�}N��uN�i]f�p}�@�Z�an�0�N�Qi?�dr��T��Lެ�U)Y�K���u��r��n����yX�Zt�_�R�b����<�u.-�Li��BG{���a�sy���ڝ��3���~��d�[�N�R�p�5������.�R���w��\t��l��r��������ia�	�����o��	C�Z������r��_�D�����SuF�J+��q}C�W�����&�����\�W@��a���{���Y���bg�]gG�J�v�f�M��90�c���_6��c4����|�����cy�tj��i����d�ζ��Η���#J���	�Ԅ���h��<A��bN�����zK���e��A4��@|��;�
����M��{�"V�$������&{�  �6���&�-�<�($��&�-�"�"��9e,�4���)��1T��1���)X�5���5Ŀ�D���)�
�:�X�9o��(?2�,���<֐�F�F�n;�$�� ��E�T����"c��8k��4T�9l;�502�f5?�hbN�u�R�lً��;d�p:���O�q���X��T�K�f���j���u���eb/�a�]�e�*������E��F$������k������ϴK��I����b��{��ҧ���h��:*�2u~�0b�$���<��9h�2���9���92���
�*��dj&�e�)��[��{��s;	�������������j;�yH�+���QS���q���bN����4���G��@�Y�UK�T�������qu��s��p��r��i��R2������-�%��*�X�V�x�`����(g��%��j4L�f�]�[����S��j�5���q�R��$��z9X�v�����j����xи���6���X�u����s�����k+�z��no��`D�P#�������L����߸���g���Z��^�����g������P�%�#�Y�8Xy�(r��P����.�?�o��s����\t�q��e��������3�~���w��w��CA���D_�R��	�z�Y.�$�W�G�3�(���V���jp��y���;��t����P�`A��z ��vw:�'"��"n��x��5��o���%�������c����wF��z��g-��E+���$I��`�j�_m�x���r\6��&W���t���D�x�5���\����p=q���~���������S������7����������e�����|��N��q���I��V���K�i�A�c!���b���W�ƻ6�ꟾ�~�w���G������\�[������d6��g���n�I��\���� ����;�D-���ʆ��p�y���A(�~o��g�����9���U�e���o���z�������������H�qg��|����hs�~���CY���=��R����l��r�!�}���nv���Z�����/*��Z���ƨ��������������*��7�a�V&-��si�h�P�"�8�v��i��c��b������_U�`���{�,��bN�g�:�|I���;��g�	��m��p�`�l�/�v$��n6��p�� %�;t��3�Y�w��jr��Q
=�K�m�Q}��9�r�,��)��7t�����=O��.�R�M��F�X����.��A��.���P4n���3���8G�'"��Kh5�k�\j�?J��>{�X�#�7o�P&��4��6ȴ�`e��Pm:�B�����Z,2�,��2�\�U���#���:T ���9?%�T��^��t���q���"��$/�������=q��"���m�����T��C��C.��.Nt�-����K����Y$��ff�^
��q���r���_�G�L]d�C���UF�X��K�Y���F�M��=q�g�P��oF��x������M��9X��R6�X�u�ě�?ҿ�#��6ȴ�=����?�T?MO�>o�{=�B^>�Z7���>����Y�?\j?� ��߀?�S�@I���1a�/�?ÿ�p?5:?�a��������?2���+?���?$�g?4��?����I�;B)�H�>�
���?	7L>o@�8���g@�ÿ�j�-B��O�?w�ٿ�f��|�����@
�X��_�?A��>ʘ/>)�>��09>��q���>aG�?<��{S?mپ�{S?���~?�V�14��B��?>��@�ݿJ�-?�-����Vׁ��?">�ܱ@ S�;Q����?��>6 ?���?΁���K�>���>vɾ�t�@�9�r��y��Þ޿�\�����]D����mS�>My�J��Zʿ�Mo= �    