"""Rebuild guide diagrams: Python + Pillow and a Japanese TrueType font.

Run from any directory. On non-Windows systems pass --font /path/to/font.
Images are static reading aids; no user data or external services are used.
"""
from pathlib import Path
import argparse
import math
from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parent


def render(name, title, labels, note, columns=3):
    rows = math.ceil(len(labels) / columns)
    width, cell_w, cell_h, gap, margin = 1140, 320, 112, 60, 30
    height = 110 + rows * 178 + 65
    im = Image.new("RGB", (width, height), "#f5f8fc")
    d = ImageDraw.Draw(im)
    d.text((margin, 24), title, font=heading, fill="#183450")
    positions = []
    for i, label in enumerate(labels):
        row, col = divmod(i, columns)
        if row % 2:
            col = columns - 1 - col
        x, y = margin + col * (cell_w + gap), 110 + row * 178
        positions.append((x, y))
        d.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=16,
                            fill="#ffffff", outline="#6790b3", width=2)
        d.ellipse((x + 13, y + 12, x + 47, y + 46), fill="#215c8e")
        d.text((x + 30, y + 28), str(i + 1), anchor="mm", font=small, fill="white")
        box = d.multiline_textbbox((0, 0), label, font=body, spacing=6)
        if box[2] - box[0] > cell_w - 62 or box[3] - box[1] > cell_h - 16:
            raise ValueError(f"Label too large in {name}: {label}")
        d.multiline_text((x + 58, y + 22), label, font=body, spacing=6, fill="#183450")
    for (x, y), (nx, ny) in zip(positions, positions[1:]):
        if y == ny:
            direction = 1 if nx > x else -1
            a = (x + cell_w + 8, y + cell_h / 2) if direction == 1 else (x - 8, y + cell_h / 2)
            b = (nx - 9, ny + cell_h / 2) if direction == 1 else (nx + cell_w + 9, ny + cell_h / 2)
        else:
            a, b = (x + cell_w / 2, y + cell_h + 8), (nx + cell_w / 2, ny - 9)
        d.line((a, b), fill="#215c8e", width=4)
        angle = math.atan2(b[1] - a[1], b[0] - a[0])
        points = [b] + [(b[0] - 13 * math.cos(angle + delta), b[1] - 13 * math.sin(angle + delta)) for delta in (-0.5, 0.5)]
        d.polygon(points, fill="#215c8e")
    d.text((margin, height - 50), note, font=small, fill="#40566d")
    im.save(OUT / name)


def roadmap():
    im = Image.new("RGB", (1140, 680), "#f5f8fc")
    d = ImageDraw.Draw(im)
    d.text((30, 24), "目的に合う手順書を選びましょう", font=heading, fill="#183450")
    d.rounded_rectangle((30, 100, 1110, 195), radius=16, fill="#215c8e")
    d.text((55, 126), "まず：操作の基本・環境構築 → 学生向け導入でFrameWorkを取得", font=body, fill="white")
    tracks = [
        (30, "計算をする方", ["外部実行環境", "準備・ビルド・実行", "ParaViewで結果を見る"]),
        (410, "プログラムを変える方", ["Gitで作業ブランチを作る", "コード拡張・テスト", "Git運用：レビュー・共有"]),
        (790, "GPEを詳しく調べる方", ["GPE / GP3D仕様書", "設定・数式・出力を調べる", "必要な章を参照する"]),
    ]
    for x, title, lines in tracks:
        d.line((x + 160, 195, x + 160, 242), fill="#215c8e", width=4)
        d.polygon([(x + 160, 252), (x + 151, 238), (x + 169, 238)], fill="#215c8e")
        d.rounded_rectangle((x, 264, x + 320, 580), radius=16, fill="white", outline="#6790b3", width=2)
        d.text((x + 17, 287), title, font=body, fill="#215c8e")
        for i, line in enumerate(lines):
            d.text((x + 17, 353 + 62 * i), line, font=small, fill="#183450")
    d.text((30, 620), "すべてを実行する必要はありません。自分の目的に合う列を選びます。", font=small, fill="#40566d")
    im.save(OUT / "reading-map.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font", default="C:/Windows/Fonts/meiryo.ttc")
    args = parser.parse_args()
    heading = ImageFont.truetype(args.font, 32)
    body = ImageFont.truetype(args.font, 23)
    small = ImageFont.truetype(args.font, 21)
    roadmap()
    render("setup-flow.png", "PCの用途を決めてから必要なソフトを準備する", [
        "PCの用途を決定\n編集・CPU・GPU・接続", "Git・エディター\n取得・編集の準備", "Python・venv\n運用と後処理の準備",
        "ビルド用の道具\nFortran・CMake等", "必要な追加部品\nMPI・FFT・CUDA", "各ソフトの版確認\n呼出し先も確認",
        "FrameWorkを取得\n自分のPCの設定へ", "最小ビルド・テスト\n基本の動作を確認", "検証ケースを実行\n対象構成の出力確認"],
        "編集・接続だけのPCは計算用の道具を省略できます。計算機ごとに準備します。")
    render("folders.png", "作業する場所とファイルの流れ", [
        "FrameWork\nソース・文書の編集元", "実行環境を生成\nResearchRuns内へ", "case.yaml\n計算条件を編集",
        "入力生成・計算\n実行環境で行う", "output\n計算結果を保存", "paraview\n表示用データを保存"],
        "矢印は作業の順序です。ソースの編集元と計算する場所を区別します。")
    render("git-changes.png", "保存した変更を共有するまで", [
        "編集・保存\nCtrl+S", "git add\n記録する変更を選ぶ", "git commit\nPCに履歴を記録",
        "git push\nGitHubへ送る", "Pull Request\n担当者に確認を依頼", "承認後にマージ\n公式mainへ取り込む"],
        "保存・commit・pushは別の操作です。Pull Requestの作成だけでは公式版は変わりません。")
    render("git-fork.png", "fork方式：自分の変更を公式へ届ける", [
        "公式リポジトリ\nupstream", "fork：GitHubに複製\n自分のorigin", "clone：PCへ取得\n作業ブランチで編集",
        "commit・push\n自分のoriginへ", "Pull Request\n公式へ確認を依頼", "承認・取り込み\n公式のmainへ"],
        "公式の更新は別途fetchしてPCのmainへ反映します。originの意味は取得方式で異なります。")
    render("run-flow.png", "計算の準備から結果確認まで", [
        "環境設計書\nモデル・保存先を決定", "dry-run\n生成計画を確認", "環境生成\n生成先へ移動",
        "case.yamlを編集\n保存して入力生成", "validate-only\n設定を検査", "build\n実行プログラムを作る",
        "run\n計算を行う", "outputを確認\nログ・時刻も確認", "後処理\n結果を図で見る"],
        "各段階で成功を確認します。エラーが出たら止まり、次の操作へ進みません。")
    render("postprocess.png", "計算結果をParaViewで表示するまで", [
        "output\nSLF ＋ meta.json", "Pythonで変換\nVTIとPVDを作る", "collection.pvd\nParaViewで開く",
        "Apply\nデータを読み込む", "密度などを選ぶ\n数値に色を付ける", "Slice\n内部の断面を見る"],
        "PVDは時刻とVTIの一覧です。別PCへ移すときはVTIを含むフォルダごとコピーします。")
    render("paraview-controls.png", "ParaViewで探す場所（画面の役割を示す模式図）", [
        "左：Pipeline Browser\n操作するデータを選ぶ", "左：Properties\n設定してApply", "中央：表示領域\n形・色を確認"],
        "上部の変数・時刻選択で表示を切り替えます。画面の配置は版によって異なります。")
    render("development.png", "プログラムを拡張するときの流れ", [
        "仕様メモ\n変更する目的を決める", "作業ブランチ\nFrameWorkで作成", "コードを変更\n設定と計算をつなぐ",
        "CMake・manifest\n必要なファイルを登録", "小さなテスト\n変更の結果を確認", "新しい実行環境\n生成・ビルド・試験",
        "文書・検証結果\n条件と結果を保存", "commit・push\n変更を共有", "Pull Request\nレビューを依頼"],
        "生成環境での試験は、必要なファイルのコピー漏れがないかも確認します。")
    render("gpe-flow.png", "GPE / GP3D：条件から結果まで", [
        "case.yaml\n人が条件を書く", "prepare\ninput.nmlを生成", "GP3D\n初期状態を作る",
        "時間発展\n方程式を解く", "SLF ＋ meta.json\n数値と格子情報を保存", "後処理・ParaView\n密度や位相を見る"],
        "普段編集する入口はcase.yamlです。input.nmlは自動生成されます。")
