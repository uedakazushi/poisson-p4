# CB generic V計算のexact QQ(V) frame

2026-11-03。CAのV1計算と分離し、CBではmodulus Vを固定しない。
既に全理想が証明されたCB source critical schemeのmonic V-axisを、

    (s0,...,s9)=(V,0,0,V/2,−V²/4,0,0,−V/2,V²/4,−V³/24)

として全てのcritical generatorに代入し、QQ[V]で零と確認する。
full ratio*f sourceを使って、その曲線の実150 tensor係数とfiber Jacobianを構成する。

QQ(V)上でGL orbitのrankは19である。
sourceのGL19と10個のfixed fiber変数を選ぶと、全29次元source chartのframeとなり、そのdeterminantは

    V^12 (V²−2)/3072

である。source-to-targetの全150×29 Jacobianはrank28、一つのkernelを厳密に保存する。
kernelの不要な共通因子を除き、次のmixed Hessian計算で次数を増やさないようにした。
これはgenericなframeの選択であり、V=0やV²=2をactual geometryの例外点と認定する結果ではない。

次のsource formulaはGL19のordered exponentialsをfull CB ratio/f/H sourceに作用させるものである。
10 fiber変数にratio方向とf方向が含まれるため、二次項

    −delta_q * sum(delta_f_i * G_i)

を全て残す必要がある。
CAのlinear fixed-fiber helperをそのままCBへ流用してはいけない。
保存したgeneric frameはこのmixed termを明記し、旧全CB physical certificateに結び付けている。

本packetはQQ(V)上の全source frame/Jacobianだけを証明する。
actual odd ideal、八regular odd constraint、cusp order、CB全完成環の主張はない。
次のbounded stepはこのframeに対するmixed kernel Hessianとactual odd firstjetである。
GL rank19なので、CAのdense single-orbit argumentをCBへ適用しない。

再現:

    .venv/bin/python code/ambient_CB_V_family_frame_20261103_v63.py --output /tmp/cb-v-frame
