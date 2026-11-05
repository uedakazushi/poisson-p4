# 六adjointのhigher-source gate：全661補正とm2の置換

2026-11-04。646 X-onlyと15 mixed d3を全て保持し、実source一次値だけを使う。
X3、m2、source/gauge高次係数は評価しない。

全661 sourceのX部分uについて、X1u=Xp0 a1を全44×8成分で確認する。
free五商での零性だけを使っていない。mixed d3はX部分零として別に保持する。
これにより独立のm2値を要求せず、18本の実pivot二次係数へ置き換えられる。

canonical m第7行はrだけの函数で、mu0=(5/2)e20、mu1=r*(保存済み44係数)である。
実cycle d1*m*Xp=0をz=0へ制限するとr*(mXp)7=0、従って(mXp)7=0。
r二乗係数により

  mu2*Xp0=−mu1*Xp1_r−(5/2)Xp2_row20,r二乗

を得る。全661でX1u=Xp0 a1なので、retained wedgeのmu2*X1u_ziはこの18 contractionと
a1_ziの積へ置換できる。m2を零にしたのではない。新しい18入力も実N2_B11/D、
既存A1/B1とnativeB11から得られる同じsource函数で、新gauge三次は不要。

constant単位blockを正規化した後、一次零sourceの二次Schur値にはD2−C1U0^-1B1だけが
現れるため、unit inverse自体の二次微分U2も要らない。一方、81pivotのfirst carry全体が
消えたとは主張しない。これを消すか、必要なadjointだけを既存cacheへbindするgateは残す。

## Native最高係数の独立感度

native rowhに現れる未知u5_even,nの係数はnに依らないので、既証の二次感度零が三次にも
同じconstant source上で適用できる。core quadraticに未知proper最高係数を入れたscalar responseは
pure tau originとnu0=0から零となる。さらにnuA,n・nuB,nを独立記号として、C4 originでの
全native140 responseを作り、同じ681 constant sourceへ運ぶ。この検査も最高係数の値ではない。
結果はreportへ保存し、零でなければそのcallbackを残す。

固定source検査は3.06秒で完了し、全661×44×8の一次補正と全780 pure-r既知係数を確認した。
nuA最高係数・nuB最高係数のnative responseは、全681 sourceでともに零だった。
u5_even,3も同じconstant係数mapが零なので、native row20最高三次inputは全て消える。
これはraw source X3を零としたのではなく、正確なconstant source functional上の相殺である。
従ってTnilから新しいnative G3/gauge3 functionalsを共同source jobへ追加する必要は無い。
旧640のlower-second積と81pivot first carryのgateは依然保持する。

旧640 row20のX3にはA0=0とBold0u=0を使えるが、A1*Bold2+A2*Bold1は保持する。
必要なBold2方向はr二乗とz0r,z4r,z5rへ限定され、全行列を要求しない。
新しい270 zz gauge cacheはexact CA-F9 base20としてpinし、N2_B11と共有する。
古い619側の五tripleと今回の三tripleを混同せず、共同source jobはまだ実行しない。
