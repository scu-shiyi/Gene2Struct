import re
from typing import Tuple, Dict, List, Set
import unicodedata


def _strip_invisibles(s: str) -> str:
    zero_width_chars = ['\u200b', '\u200c', '\u200d', '\ufeff']
    for char in zero_width_chars:
        s = s.replace(char, '')
    return s


def _sanitize_id_strict(s: str) -> str:
    s = unicodedata.normalize('NFC', s)
    s = _strip_invisibles(s)
    s = s.replace('\u00A0', ' ')  # NBSP -> space
    cleaned_id = re.sub(r'[\s\.]+', '_', s, flags=re.UNICODE)
    cleaned_id = re.sub(r'[^\w]', '_', cleaned_id)
    cleaned_id = re.sub(r'_{2,}', '_', cleaned_id)
    cleaned_id = cleaned_id.strip('_')
    return cleaned_id

def parse_tree(tree_file):
    translate_dict: Dict[str, str] = {}
    species_data: List[Tuple[str, str]] = []  # 存储 (ID, Name)
    species_map: Dict[str, str] = {}  # 用于 Newick 格式生成 species_data
    species_id_counter: int = 1
    in_translate_block: bool = False
    has_translate: bool = False

    try:
        with open(tree_file, "r", encoding="utf-8") as file:
            lines = [line.strip() for line in file.readlines()]
            for line in lines:
                if not line:
                    continue

                # --- 1. 识别 TRANSLATE 块 ---
                if "TRANSLATE" in line.upper():
                    has_translate = True
                    in_translate_block = True
                    continue

                # --- 2. 解析 TRANSLATE 块内容 ---
                if in_translate_block:
                    if ";" in line:
                        in_translate_block = False
                        continue
                    if line:
                        # 假设格式为 "ID Name," 或 "ID 'Name',"
                        parts = line.replace(";", "").split()
                        if len(parts) >= 2:
                            species_id, species_name = parts[0], parts[1].replace("'", "").strip(",")
                            cleaned_name = _sanitize_id_strict(species_name)
                            translate_dict[species_id] = cleaned_name
                            species_data.append((species_id, species_name))
                        continue

                # --- 3. 识别和处理 Newick/Nexus 树字符串 ---
                # 典型的 Newick 格式特征：包含 (), :, 和 ;
                if "(" in line and ")" in line and ":" in line and ";" in line:
                    # 清理 ete3/FigTree 附带的额外信息
                    line = re.sub(r"\[&&NHX:[^\]]+\]", "", line)
                    line = re.sub(r"\[&[^\]]*\]", "", line)
                    # 移除 Nexus 格式中的 'tree NAME = ' 前缀
                    current_tree = re.sub(r"^\s*tree\s+\w+\s*=\s*", "", line, flags=re.IGNORECASE).strip()

                    if has_translate:  # Nexus 格式，需要翻译

                        # 替换 Newick 字符串中的 ID 为实际物种名
                        temp_tree = current_tree
                        for id, name in translate_dict.items():
                            # 使用 r'\b' (单词边界) 确保只替换完整的 ID
                            temp_tree = re.sub(rf'\b{id}\b', name, temp_tree)
                        newick_tree = temp_tree

                    else:  # 纯 Newick 格式，直接提取物种名
                        newick_tree = current_tree

                        # 提取物种名 (可能包含 ID, 名称或别名，只要后面跟着 :距离)
                        # 正则表达式：匹配 (字母数字、点、下划线、横杠组成的) 物种名, 后面跟着 :距离
                        matches = re.findall(r"([A-Za-z][A-Za-z0-9._-]*):[\d.]+", newick_tree)

                        # 为每个物种生成 ID (如果需要 species_data)
                        for species_name in matches:
                            cleaned_name = _sanitize_id_strict(species_name)
                            if cleaned_name and cleaned_name not in species_map:
                                species_map[cleaned_name] = str(species_id_counter)
                                species_data.append((str(species_id_counter), cleaned_name))
                                species_id_counter += 1

                    break

            # 从 (ID, Name) 列表中提取所有唯一的物种名称
            species_set = set([i[1] for i in species_data])
            return species_set, newick_tree

    except Exception as e:
        raise ValueError(f"Failed to read tree file: {e}")

if __name__ == "__main__":
    species_set, newick_tree = parse_tree('//evolution/out/rbcl.treefile')
    print([i for i in species_set])
    print(newick_tree)