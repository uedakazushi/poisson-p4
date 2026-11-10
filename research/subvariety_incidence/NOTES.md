# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v6。数学的証明の全文はpapers/poisson_subvariety_incidence_ja.texと共有papers/quadric_log_orthogonality.tex。再開時はRESEARCH_POLICY.mdを最初に読む。

## 非退化三次元qUの共鳴Log全層別

V=U1⊕L⊕W、ξ=(1,1,-2)としたPoisson源 P0=2ξ∧R+ξ(qU)J。qU非退化の場合、H=qU|U1は階数2か1（階数0ならqU退化）。

- rankH=2: qU=H(u)+2zb(u)+cz²、b≠0ならstrict stabilizerはCRだけで全延長A1（rank2）、b=0ならCR⊕so(H)でJacobi理想(v(t+2))、rank4因子Q+V(z)+V(cz²+w0²+w1²)。c=0でもqU非退化ならb非等方となり、追加対称性はない。
- rankH=1: 重みと可換なGL(U1)でqU=u0²+2u1z+cz²へ正規化できる。c≠0ならstabilizer CRで全延長A1、c=0なら追加場A=-u1∂u1+z∂zが現れる。KA=-q∂u1∧∂z+E∧A、全Jacobi idealは(v(t-1))、rank4因子はQ+V(u0)+V(u0²+w0²+w1²)。他の未確認な厚みはない。
- 直交性wall b=0 と階数1の等方性wall c=0 はこのmarked二次形式族でいずれも余次元2。個別の147成分同一視や非包含をこれだけで主張しない。

厳密なrank1算術はcode/subvariety_incidence/check_log_rankone.py（c=0,1でstrict so5 rank8/9、Jacobi GB=(v(t-1))、Pfaffian）。一般定理の証明は零スキームZ(P0)の線形包がWとUを回復すること、so(qU)とξのcentralizerの線形代数、qU正規形にある。

## 残課題と方向転換条件

qU退化の共鳴Log境界、special Affとの交叉とrank4枝の147成分への帰属は未解決。これらが不必要に細かい分類へ向かうなら、4次del Pezzoの延長例外、P3の横断的Poisson対称性、147成分の共通幾何へ方向を変更する。元の一般Q3 Rat/Log/Affの延長、smooth cubic一意剛性、一般del Pezzo剛性は現行原稿に保存。新しい主張の文献上の優先権は未確認。

保存・再現：正本GitHub main、TeX＋7 producer/7 JSON、PDF非追跡、GitHub Actionsなし、Libraryはpush失敗時だけ。
