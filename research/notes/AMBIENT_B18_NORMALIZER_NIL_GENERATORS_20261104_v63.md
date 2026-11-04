# normalized E4 の全 nilideal：明示二生成と完全な module 関係

2026-11-04。対象は実正規化比較環

    E4=S9/(dJ),
    S9=k[[u1,...,u5,a,d,x,t]],
    J=(u1,...,u5,Q0,Q1,(a+x)*(x-a-d)^2)。

Q0,Q1 は既存 actual source の11項・27項多項式であり、full core CI と五直線
支持は独立監査済みである。元の D0 の nilideal やその生成数は主張しない。
特に nil-cotangent の次元から principal ideal と推測する方法は使わない。

## 1. 全 minimal primes から radical を計算する

受理済みの五 core line のうち、N,P は d=0に入り、C,R,L は入らない。
従って E4 の minimal prime は (d) と (u5)+I_C、(u5)+I_R、(u5)+I_L。
embedded N/P を minimal intersection に重複させない。これらを intersect すれば

    sqrt(dJ)=(d) intersect (u5+I3),
    I3=I_C intersect I_R intersect I_L。

以下の有理 GL4 座標を使う。

    c=(d+x)/2, r=(d-x-2a)/6, l=(d-x+a)/3,
    h=a-d+x-2t。

逆変換は a=-2r+l、d=c+r+l、x=c-r-l、t=(-4r-l-h)/2。
実三直線 C,R,L はそれぞれ c,r,l 軸となり、いずれも h=0。
従って

    I3=(h,c*r,c*l,r*l)。

この quotient は被約三軸環であり、d=c+r+l はどの minimal prime にも入らない。
ゆえに d は S9/(u5+I3) の非零因子であり

    sqrt(dJ)=d*(u5+I3)。

したがって全 nilideal は正確に

    Nil(E4)=d*(u5+I3)/d*(u5+Jcore),
    Jcore=(Q0,Q1,(a+x)*(x-a-d)^2)。

S9 自体では d が非零因子なので、d を掛ける前の module と同型であり、

    Nil(E4) isomorphic to I3/Jcore

となる。この module は u5 に annihilate されるが、E4 自体の u5 を除いていない。
全五つの d*u_j 関係は環の定義に保持する。

## 2. 実二次式で二つの quad generator を消す

g=(h,c*r,c*l,r*l) とし、Q0,Q1,f=(a+x)*(x-a-d)^2 のそれぞれを g の結合へ
固定 monomial division で書く。その全係数を保存する。I3/m*I3 は四次元で、
Q0,Q1 の class は三つの quad class の中で rank2、f の class は零となる。
source 多項式を実際に代入した厳密 QQ 行列は、行を(h,cr,cl,rl)、列を(Q0,Q1,f)
とすると

    [[0,0,0],[32,0,0],[11,0,0],[-9,-90,0]]

である。cr と rl の2x2 minor は -2880、残す quad は q=c*l。

origin で単元となる2x2 quad-coefficient minor Delta を選び、その adjugate
だけで二つの quad を Q0,Q1、h、残る quad q の結合へ解く。
各項には同じ既知単位 Delta の分母だけを用いる。全四つの恒等式

    Delta*g = E*(h,q) + T*(Q0,Q1)

を literal polynomial identity として保存・検査する。よって全 nilideal は
d*h、d*q の二つで生成される。二つの独立な class が I3/(Jcore+m*I3) に
残るので、この生成数2は最小である。実 original core 座標で書けば

    n1=d*(a-d+x-2t),
    n2=d*(d+x)*(a+d-x)/6。

新しい分母は全体で単位 Delta のみであり、Delta(0)=-2880を全係数で確認した。
GL4変換の行列式は1/3。正規化比較環における結果であって
元 D0 の nil generator 数ではない。

## 3. 完全 module presentation

I3 の全第一 syzygy は次の五本で生成される。

    (c*r)*e_h-h*e_cr,
    (c*l)*e_h-h*e_cl,
    (r*l)*e_h-h*e_rl,
    l*e_cr-r*e_cl,
    l*e_cr-c*e_rl。

完全性の証明：hを零にすると三つの cross quad の syzygy は後二本の
Hilbert–Burch columns で生成される。それらを引けば各 quad 成分が h の
倍数となり、h を正則環内で除いて前三本の crossing relations を得る。
GB や有限 jet の閉包推測を必要としない。

I3/Jcore はこの五本と Q0,Q1,f の三本を加えた4-generator/8-relation module。
前節の単位2x2消去で Q0,Q1 を消すと、(h,q) の **完全2x6 presentation** を
core S4=k[[c,r,l,h]] 上で得る。六本それぞれについて、その像が Q0,Q1,f の
どの結合に等しいかを全 QQ polynomial cofactor まで保存する。

係数環を full S9 または E4 とすると u5 の作用も必要である。従って上の
六本に u_j*e_k（j=1,...,5、k=1,2）の十本を付けた **完全2x16 presentation**
を使う。2x6だけを full S9/E4 の全関係としない。六本/十六本の最小性は
ここでは主張しない。明示した生成元と S9 作用が nilideal の積も一意に指定する。

計算は実 source 三多項式と小さい有理線形代数・polynomial 恒等式だけである。
新しい source 点・jet・一般体・GB を使わない。保存された新多項式・matrix係数は
非零827項、一回の計算は0.95秒以内で完了した。具体的な q、Delta、二 nil
生成元、全変換・cofactor・六関係は同名 report と証明書に保存する。

二生成元とも m_E4^2 に入るので、ここでの nil-cotangent image は零である。
それでも nilideal の最小生成数は2である。この差からも、別の環 D0 の
nil-cotangent image の次元だけで principal nilideal と推論してはいけない。
