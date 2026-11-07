# μ五 marked H+quartic 四候補の全周囲成分性の独立確認

2026-11-07、PRIMARY E。別実装の厳密Fractionで四件のfullambient componenthoodを確認した。
全sourceの一23軌道還元・全marked geometry・既111及びμ六九候補との相違は別作者監査を要する。
A5 smooth hyperflexの二normal sourceには高次計算を加えていない。

`classification_first_essential_quintic_review.py --mu5-marked-ambient` は、
μ六で用いた全係数アルゴリズムを共通packet engineへ追記したmodeである。
旧十一FunctionDefのASTは全て不変。旧a818c890 source全bytes及び実行jobは
`inputs/marked_quartic_reviewer_a818c890_frozen_20261107.json` に保存した。
producer・legacy・CASはimportしない。契約 `jobs/plane_hyperplane_quartic_mu5_ambient_review.json`
のgate後、一回実行し8.747秒で全PASS。

A4+A1 node inflexion、A4+A1 smooth hyperflex、D5 smooth hyperflex、2A2+A1 smooth hyperflexについて、
全六磁気sourceから同一canonical点を復元し、全350Jacobi、全350×150 Dとrank101、
Euler/GL/scalar/磁気の全57列とrank48、全Dとの積零、非零minor、全126実D及びessential rank三を確認した。
実48列＋一本normal＋101IFT座標の全150直和detを独立に再構成した。
全350二次Jacobi、全150 D列を消すfunctional及び非零pairingも一致する。

真の102次元cutでIFT後の唯一normalが非零二次leadingを持つため被約支持は原点であり、
実48方向の下界と合わせraw局所次元48、projective23軌道閉包が全成分となる。
smoothnessや非被約厚さ、高次局所環の主張はしない。
直和の列操作による49行縮小とKrullの論理は
`PLANE_HYPERPLANE_QUARTIC_MU6_AMBIENT_SEPARATE_REVIEW.md` と同一である。

report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_AMBIENT_REVIEW_20261107.json` は算術証拠であり、
全curve/source分類や既成分との相違を算術だけから結論しない。
それらの独立確認がそろって初めて四件を採録する。全G/Eと上界は依然未証明である。
