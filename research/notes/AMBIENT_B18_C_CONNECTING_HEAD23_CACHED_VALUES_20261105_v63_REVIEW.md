# head23 の一 cache 合成と十四出力：独立 VALUE 検証

独立 VALUE 検証は PASS した。報告 SHA256 は
`2b2669c98b68953eb6e2826db62012a4e20b1f1ceffb19ab01fd00d059f41689`。
原 VALUE は `7d51e4a1…`。原2,194 record、同じ supplied105×10 branch、
五つの root 低次係数、単位逆の左右積、全十四出力と全項別積が一致した。

producer を import せず、与えられた r,s の全十係数を binomial 公式で
s=r+ell/5 に移した。環は QQ[r,ell]/(r²,ell³) で、係数は単項式対を key とする
辞書で保持した。唯一の単位逆は、一次以上の nilpotent 部分 t を使って
(1+t)⁻¹=1−t+t²−t³ と再構成した。t⁴=0 と左右逆の両方を確認した。
これは原 worker の係数順の三角 recurrence とは異なる計算法である。

この逆計算法の追加事前確認として、単位入力の517-node閉包だけに正値 majorant を伝播した。
原 root の signed 正次数値を求める前に、geometric powers と左右逆診断が同じ41 bit上界に
入ることを示し、root83 bit・post207 bitの演算前 guard と全所有量を保存した。
2,194 record の実 cache は last-use120 nodeを保持し、source gate と一致した。
provided105 の変換 fingerprint も原記録と同一だった。

その後、十四多項式を項順・積順とも逆向きに評価した。
全ての項別積を元の記録と照合し、出力だけの一致に依存していない。
得た一つの新しい既存 cache 合成係数は

  bB63[r ell²] = 2039/360

である。head23 の canonical k2 座標は 1331/9216、h0 の低 ellH 座標は2/25。
従って指定された単純 target minor は

  (2/25)(1331/9216) = 1331/115200 ≠ 0

である。

受理済みの全射 M7=J2²⊕k³→G4=J2⊕k² を使う。
この minor が非零なので、Soc M の像は rG と二単純 target を含み三次元以上である。
三次元 kernel が半単純ならその像は二次元となって矛盾する。
従って kernel の r-rank は1で、実型は J2(r)⊕k となる。
これは honest M7 の kernel であり、Kbar/IKbar ではない。
他の head の高座標、元 K の torsion 指数、元 H0、nil 三乗への転用は行わない。

実行は0.063秒、所有上界37,142、最大演算前114 bitであった。
独立 job も source83→207、hard2048、60秒、2 GiB、60,000 slots、RAW2 MBの範囲を守った。
新しい物理・IFT・gauge・jet solve、他の source query はない。

再現コードは `code/ambient_B18_C_connecting_head23_values_review_20261105_v63.py`。
数学 SOURCE と cache/code SOURCE の両 peer を入力として pin した。
凍結出力先への再実行は拒否するため、別作業コピーで再検証すること。
