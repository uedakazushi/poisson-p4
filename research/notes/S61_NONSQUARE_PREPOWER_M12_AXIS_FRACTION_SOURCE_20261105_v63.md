# 追加axis source評価のFraction temporary上界

2026-11-05。新しい符号付きaxis値は計算しない。
保存source profileの全214中間、追加列/行、row LCD、整数row L1の最大状態は118bit。
各状態の既約numeratorはそのL1上界以下、denominatorは記録した共通Dの約数である。
Fractionの一回の加減算を未約分cross積で行っても、numeratorは二つの118bit積の和、
denominatorは118bit積で足りる。poly積再生のpartial remainderのL1は二倍以下なので、
最大temporaryは2*118+2=238bitとなる。全入出力のLCD/a shift/QQ scaleも同じ範囲に入る。
よってsource axis prepにも2048bitを越える仮定は不要。

全q214値を解放せず保持してもaxis支持の和は1309。key/numerator/denominatorで3927slot。
source metadata、旧579、選択decodeの一cell、二十三追加poly、native整数input、
callback bufferを別に数えた全prep予算は二万slot以下で、既存33275に収まる。
実装時には全parameterの巨大JSONをdecodeせず、必要subtree/cellだけを読む。
新axis値・determinant・rankは未実行。
