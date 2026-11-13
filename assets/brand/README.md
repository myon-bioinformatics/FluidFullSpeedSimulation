# FluidFullSpeedSimulation のブランド画像

提供シート（1536×1024 PNG、sha256 `88333fab2f059083533d1cd8f2eb130b3555defe74a6961597ffadb67c6f2cf2`）から、日本語キャプションを含めずに 3 点を切り出した。原本シートはリポジトリに置かない。座標はシート左上原点・端点含む。

| ファイル | 役割 | 大きさ | 切り出し（左, 上）–（右, 下） | sha256 |
|---|---|---|---|---|
| [icon.png](icon.png) | アプリアイコン。README 上部に縮小表示する | 512×512 | (17, 104)–(430, 523)（414×420 を左右 3 px の `#FBFAF9` で 420×420 に揃え、LANCZOS で拡大） | `46355eaa9580ac8d085c68d5d39610189b2be82b4699d519dded271d16cf3450` |
| [favicon.png](favicon.png) | 専用ファビコン（竹と水）。利用ガイド（GitHub Pages）のタブアイコン | 64×64 | (148, 663)–(300, 814)（153×152 を 153×153 に揃えてから LANCZOS で縮小。アイコンの縮小ではない） | `bd7a8c88a7e4a8557005e471d6ff24407ed96e67ddb61f6d24fe917a79646668` |
| [favicon-32.png](favicon-32.png) | 同上 32×32 | 32×32 | ファビコンと同じ切り出し | `d341f380622c41e8205a128cbed95ceadc448649f53267a28ede07bcfe135f01` |
| [favicon-16.png](favicon-16.png) | 同上 16×16 | 16×16 | ファビコンと同じ切り出し | `f824be863192df3cae0eb422307f0d0baa3d67d0ca51c04a1f481964e19e1a0d` |
| [favicon.ico](favicon.ico) | 同上 ICO（16 / 32 / 64 を格納）。Pages ではサイト根の `/favicon.ico` にも複製する | ICO | ファビコンと同じ切り出し | `2d74bca3c1e68964efe0907649ef58859196331506c2c57e481b4b424a5d46eb` |
| [concept-diagram.png](concept-diagram.png) | 概念図パネル。README の概要でクリック拡大する | 1076×1010 | (454, 6)–(1529, 1015)（パネル枠 (456, 8)–(1527, 1013) に 2 px の余白） | `b8d0517a027004d0e97e30c5fd28b2ec139e3a32ff44d7140c588694bd52b3a6` |

利用ガイドの `docs/guide.page.json` には `favicon` だけを渡す。`icon` を渡すとページ本文先頭に alt のない 80×80 画像が出るため、入れていない。
