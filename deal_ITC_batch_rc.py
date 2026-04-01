##Author：cheng hao
##Email: chenghao8425@163.com
##微信公众号：材料计算学习笔记
from pathlib import Path
from typing import Dict, Any, List, Tuple
from contextlib import redirect_stdout
from datetime import datetime
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


# =========================
# 1) 从温度原始文件提取平均温度分布
# =========================
def process_temperature_data(
    input_file,
    output_temperature_file,
    output_mean_file,
    figure_file,
    num_blocks=50,
    num_chunks=50,
    position_factor=200.0,
):
    with open(input_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    temperature_data = []
    chunk_positions = []
    current_block_temperatures = []

    for line in lines:
        parts = line.strip().split()

        # LAMMPS chunk 输出常见格式：3列行为 block 分隔信息，4列行为 chunk 数据
        if len(parts) == 3 and parts[0].isdigit():
            if current_block_temperatures:
                temperature_data.append(current_block_temperatures)
                current_block_temperatures = []
            continue

        if len(parts) == 4 and parts[0].isdigit():
            if len(chunk_positions) < num_chunks:
                chunk_positions.append(float(parts[1]) * position_factor)
            current_block_temperatures.append(float(parts[3]))

    if current_block_temperatures:
        temperature_data.append(current_block_temperatures)

    if not temperature_data:
        raise ValueError(f"未从温度文件中读取到有效数据: {input_file}")

    temperature_df = pd.DataFrame(temperature_data[:num_blocks]).T
    temperature_df.insert(0, "Position", chunk_positions)
    temperature_df["Position"] = temperature_df["Position"].round().astype(int)

    Path(output_temperature_file).parent.mkdir(parents=True, exist_ok=True)
    temperature_df.to_csv(output_temperature_file, sep='\t', index=False, header=False)

    mean_temperatures = temperature_df.iloc[:, 1:].replace(0, pd.NA).mean(axis=1)
    mean_df = pd.DataFrame({"Position": temperature_df["Position"], "Mean_Temperature": mean_temperatures})
    mean_df.to_csv(output_mean_file, sep='\t', index=False, header=False)

    plt.figure(figsize=(8, 6))
    plt.scatter(mean_df["Position"], mean_df["Mean_Temperature"], alpha=0.8, label='Temperature Data')
    plt.xlabel('Position (Å)')
    plt.ylabel('Mean Temperature (K)')
    plt.title('Mean Temperature vs Position')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_file, dpi=300)
    plt.close()

    return mean_df


# ===================================
# 2) 热流线性拟合 -> 平均热功率（W）
# ===================================
def compute_heat_power(
    heat_flux_file,
    figure_file,
    time_ps_factor=0.0005,
):
    data = np.loadtxt(heat_flux_file, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 3:
        raise ValueError(f"热流文件列数不足，至少需要3列: {heat_flux_file}")

    time_ps = data[:, 0] * time_ps_factor
    f_hot = data[:, 1]
    f_cold = data[:, 2]

    slope_hot, intercept_hot, _, _, _ = linregress(time_ps, f_hot)
    slope_cold, intercept_cold, _, _, _ = linregress(time_ps, f_cold)

    factor = 1.60218e-19 / 1.0e-12  # eV/ps -> W
    slope_hot_w = slope_hot * factor
    slope_cold_w = slope_cold * factor
    P_avg = (abs(slope_hot_w) + abs(slope_cold_w)) / 2.0

    plt.figure(figsize=(10, 6))
    plt.plot(time_ps, f_hot, 'r.', alpha=0.6, label='Hot region data')
    plt.plot(time_ps, f_cold, 'b.', alpha=0.6, label='Cold region data')
    plt.plot(time_ps, slope_hot * time_ps + intercept_hot, 'r-', label=f'Hot fit: {slope_hot:.4f} eV/ps')
    plt.plot(time_ps, slope_cold * time_ps + intercept_cold, 'b-', label=f'Cold fit: {slope_cold:.4f} eV/ps')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy Exchange (eV/ps)')
    plt.title('Energy Exchange Linear Fit')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_file, dpi=300)
    plt.close()

    return {
        'P_avg_W': P_avg,
        'slope_hot_eV_ps': slope_hot,
        'slope_cold_eV_ps': slope_cold,
        'slope_hot_W': slope_hot_w,
        'slope_cold_W': slope_cold_w,
        'intercept_hot': intercept_hot,
        'intercept_cold': intercept_cold,
    }


# ===================================
# 3) 双区间温度线性拟合 -> 界面温差 dT
# ===================================
def fit_temperature_difference(
    mean_temp_file,
    range1: Tuple[float, float],
    range2: Tuple[float, float],
    mid_x: float,
    figure_file,
):
    data = pd.read_csv(mean_temp_file, sep='\t', header=None, names=['x', 'T']).dropna()
    data = data[data['T'] > 0]
    if len(data) < 4:
        raise ValueError(f"平均温度数据不足: {mean_temp_file}")

    m1 = (data['x'] >= range1[0]) & (data['x'] <= range1[1])
    m2 = (data['x'] >= range2[0]) & (data['x'] <= range2[1])
    if m1.sum() < 2:
        raise ValueError(f"range1 区间内数据点不足: {range1}")
    if m2.sum() < 2:
        raise ValueError(f"range2 区间内数据点不足: {range2}")

    x1 = data.loc[m1, 'x']
    T1 = data.loc[m1, 'T']
    x2 = data.loc[m2, 'x']
    T2 = data.loc[m2, 'T']

    s1, i1, _, _, _ = linregress(x1, T1)
    s2, i2, _, _, _ = linregress(x2, T2)

    T_left_mid = s1 * mid_x + i1
    T_right_mid = s2 * mid_x + i2
    dT = abs(T_left_mid - T_right_mid)

    plt.figure(figsize=(10, 6))
    plt.scatter(data['x'], data['T'], color='blue', alpha=0.7, label='Temperature Data')
    plt.plot(x1, s1 * x1 + i1, color='red', label=f'Fit range 1: slope={s1:.4f} K/Å')
    plt.plot(x2, s2 * x2 + i2, color='green', label=f'Fit range 2: slope={s2:.4f} K/Å')
    plt.axvline(mid_x, color='purple', linestyle='--', label=f'Interface x={mid_x}')
    plt.xlabel('Position (Å)')
    plt.ylabel('Mean Temperature (K)')
    plt.title('Temperature Linear Fit for ITC')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(figure_file, dpi=300)
    plt.close()

    return {
        'deltaT_K': dT,
        'range1_slope_K_per_A': s1,
        'range1_intercept': i1,
        'range2_slope_K_per_A': s2,
        'range2_intercept': i2,
        'T_left_at_mid': T_left_mid,
        'T_right_at_mid': T_right_mid,
    }


# =========================
# rc.in 配置读取
# =========================
def parse_scalar(text: str):
    s = text.strip()
    low = s.lower()
    if low in {'true', 'yes', 'on'}:
        return True
    if low in {'false', 'no', 'off'}:
        return False
    if low in {'none', 'null'}:
        return None

    try:
        if s.startswith('0') and len(s) > 1 and s[1].isdigit() and not s.startswith('0.'):
            raise ValueError
        return int(s)
    except ValueError:
        pass

    try:
        return float(s)
    except ValueError:
        return s


def parse_value(text: str):
    s = text.strip()
    if ',' in s:
        return [parse_scalar(part) for part in s.split(',') if part.strip()]
    return parse_scalar(s)


def load_rc_config(config_path: str) -> Dict[str, Any]:
    cfg: Dict[str, Any] = {
        'defaults': {},
        'case_overrides': {},
    }
    current_section = 'global'
    current_case = None

    with open(config_path, 'r', encoding='utf-8') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue

            if line.startswith('[') and line.endswith(']'):
                section = line[1:-1].strip()
                current_case = None
                if section.lower() == 'defaults':
                    current_section = 'defaults'
                elif section.lower().startswith('case:'):
                    current_section = 'case'
                    current_case = section.split(':', 1)[1].strip()
                    cfg['case_overrides'].setdefault(current_case, {})
                else:
                    current_section = 'global'
                continue

            if '=' not in line:
                raise ValueError(f"配置文件格式错误，缺少 '=' : {raw_line.strip()}")

            key, value = line.split('=', 1)
            key = key.strip()
            value = parse_value(value)

            if current_section == 'global':
                cfg[key] = value
            elif current_section == 'defaults':
                cfg['defaults'][key] = value
            elif current_section == 'case':
                cfg['case_overrides'][current_case][key] = value

    cfg.setdefault('scan_mode', 'list')
    cfg.setdefault('case_dirs', [])
    cfg.setdefault('batch_summary_name', 'batch_itc_results_summary.csv')

    defaults = cfg['defaults']
    defaults.setdefault('temp_input_name', 'temp.CBNT.txt')
    defaults.setdefault('heat_flux_input_name', 'heat_flux_exchange.txt')
    defaults.setdefault('extracted_temp_name', 'extracted_temperature_data.txt')
    defaults.setdefault('mean_temp_name', 'mean_temperature_data.txt')
    defaults.setdefault('mean_temp_figure_name', 'mean_temperature_vs_position.png')
    defaults.setdefault('heat_flux_figure_name', 'energy_flux_fit.png')
    defaults.setdefault('temp_fit_figure_name', 'temperature_fit.png')
    defaults.setdefault('summary_name', 'result_summary.csv')
    defaults.setdefault('log_name', 'ITC_results_log.txt')
    defaults.setdefault('time_ps_factor', 0.0005)

    required_top = ['root_dir']
    required_defaults = [
        'num_blocks', 'num_chunks', 'position_factor',
        'fit_range1', 'fit_range2', 'mid_x'
    ]

    for key in required_top:
        if key not in cfg:
            raise ValueError(f"配置文件缺少顶层参数: {key}")
    for key in required_defaults:
        if key not in defaults:
            raise ValueError(f"配置文件 [defaults] 缺少参数: {key}")

    for key in ['fit_range1', 'fit_range2']:
        if not isinstance(defaults[key], list) or len(defaults[key]) != 2:
            raise ValueError(f"{key} 必须写成两个值，例如: {key} = 35, 92")

    if 'area_m2' not in defaults and 'area_A2' not in defaults:
        raise ValueError("配置文件 [defaults] 需要提供 area_A2 或 area_m2 之一。")

    return cfg


# =========================
# 其他工具
# =========================
def merge_case_config(defaults: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged = defaults.copy()
    merged.update(override or {})
    return merged


def discover_case_dirs(root_dir: Path, mode: str, listed_cases: List[str]) -> List[Path]:
    if mode == 'list':
        return [root_dir / str(name) for name in listed_cases]
    if mode == 'auto':
        return sorted([p for p in root_dir.iterdir() if p.is_dir()])
    raise ValueError("scan_mode 只能是 'list' 或 'auto'")


def resolve_area_m2(cfg: Dict[str, Any]) -> float:
    if 'area_m2' in cfg and cfg['area_m2'] is not None:
        return float(cfg['area_m2'])
    return float(cfg['area_A2']) * 1.0e-20


# =========================
# 单个 case 处理
# =========================
def process_one_case(case_dir: Path, cfg: Dict[str, Any]) -> Dict[str, Any]:
    temp_input = case_dir / cfg['temp_input_name']
    heat_flux_input = case_dir / cfg['heat_flux_input_name']

    extracted_temp_file = case_dir / cfg['extracted_temp_name']
    mean_temp_file = case_dir / cfg['mean_temp_name']
    mean_temp_fig = case_dir / cfg['mean_temp_figure_name']
    heat_flux_fig = case_dir / cfg['heat_flux_figure_name']
    temp_fit_fig = case_dir / cfg['temp_fit_figure_name']
    summary_file = case_dir / cfg['summary_name']
    log_file = case_dir / cfg['log_name']

    if not temp_input.exists():
        raise FileNotFoundError(f"缺少温度文件: {temp_input}")
    if not heat_flux_input.exists():
        raise FileNotFoundError(f"缺少热流文件: {heat_flux_input}")

    with open(log_file, 'w', encoding='utf-8') as f:
        tee = Tee(sys.stdout, f)
        with redirect_stdout(tee):
            print("==== Interface Thermal Conductance (ITC) Calculation Log ====")
            print(f"Timestamp: {datetime.now()}")
            print(f"Case dir: {case_dir}\n")

            process_temperature_data(
                input_file=temp_input,
                output_temperature_file=extracted_temp_file,
                output_mean_file=mean_temp_file,
                figure_file=mean_temp_fig,
                num_blocks=cfg['num_blocks'],
                num_chunks=cfg['num_chunks'],
                position_factor=cfg['position_factor'],
            )
            print(f"[1] Temperature profile generated: {mean_temp_file}")

            heat_result = compute_heat_power(
                heat_flux_file=heat_flux_input,
                figure_file=heat_flux_fig,
                time_ps_factor=cfg['time_ps_factor'],
            )
            print(f"[2] slope_hot = {heat_result['slope_hot_eV_ps']:.6f} eV/ps ({heat_result['slope_hot_W']:.6e} W)")
            print(f"[2] slope_cold = {heat_result['slope_cold_eV_ps']:.6f} eV/ps ({heat_result['slope_cold_W']:.6e} W)")
            print(f"[2] P_avg = {heat_result['P_avg_W']:.6e} W")

            temp_result = fit_temperature_difference(
                mean_temp_file=mean_temp_file,
                range1=tuple(cfg['fit_range1']),
                range2=tuple(cfg['fit_range2']),
                mid_x=cfg['mid_x'],
                figure_file=temp_fit_fig,
            )
            print(f"[3] range1 slope = {temp_result['range1_slope_K_per_A']:.6f} K/Å")
            print(f"[3] range2 slope = {temp_result['range2_slope_K_per_A']:.6f} K/Å")
            print(f"[3] deltaT at x={cfg['mid_x']} Å = {temp_result['deltaT_K']:.6f} K")

            area_m2 = resolve_area_m2(cfg)
            itc = heat_result['P_avg_W'] / (temp_result['deltaT_K'] * area_m2)
            print(f"[4] Area = {area_m2:.6e} m^2")
            print(f"[4] ITC = {itc:.6e} W/m^2·K ({itc*1e-6:.6e} MW/m^2·K)")

    area_m2 = resolve_area_m2(cfg)
    result = {
        'case': case_dir.name,
        'status': 'ok',
        'P_avg_W': heat_result['P_avg_W'],
        'deltaT_K': temp_result['deltaT_K'],
        'ITC_W_m2K': itc,
        'ITC_MW_m2K': itc * 1e-6,
        'slope_hot_eV_ps': heat_result['slope_hot_eV_ps'],
        'slope_cold_eV_ps': heat_result['slope_cold_eV_ps'],
        'slope_hot_W': heat_result['slope_hot_W'],
        'slope_cold_W': heat_result['slope_cold_W'],
        'intercept_hot': heat_result['intercept_hot'],
        'intercept_cold': heat_result['intercept_cold'],
        'range1_slope_K_per_A': temp_result['range1_slope_K_per_A'],
        'range1_intercept': temp_result['range1_intercept'],
        'range2_slope_K_per_A': temp_result['range2_slope_K_per_A'],
        'range2_intercept': temp_result['range2_intercept'],
        'T_left_at_mid_K': temp_result['T_left_at_mid'],
        'T_right_at_mid_K': temp_result['T_right_at_mid'],
        'fit_range1_start': cfg['fit_range1'][0],
        'fit_range1_end': cfg['fit_range1'][1],
        'fit_range2_start': cfg['fit_range2'][0],
        'fit_range2_end': cfg['fit_range2'][1],
        'mid_x_A': cfg['mid_x'],
        'area_A2': cfg.get('area_A2', area_m2 * 1e20),
        'area_m2': area_m2,
        'position_factor': cfg['position_factor'],
        'num_blocks': cfg['num_blocks'],
        'num_chunks': cfg['num_chunks'],
        'time_ps_factor': cfg['time_ps_factor'],
    }

    pd.DataFrame([result]).to_csv(summary_file, index=False)
    return result


# =========================
# 主程序
# =========================
def main(config_path='rc_itc.in'):
    cfg = load_rc_config(config_path)
    root_dir = Path(str(cfg['root_dir'])).resolve()
    root_dir.mkdir(parents=True, exist_ok=True)

    defaults = cfg['defaults']
    scan_mode = str(cfg.get('scan_mode', 'list'))
    listed_cases = cfg.get('case_dirs', [])
    if isinstance(listed_cases, str):
        listed_cases = [listed_cases]
    case_overrides = cfg.get('case_overrides', {})
    summary_csv_name = str(cfg.get('batch_summary_name', 'batch_itc_results_summary.csv'))

    case_dirs = discover_case_dirs(root_dir, scan_mode, listed_cases)
    if not case_dirs:
        raise ValueError(f"在 {root_dir} 下没有找到待处理文件夹。")

    all_results = []
    for case_dir in case_dirs:
        case_cfg = merge_case_config(defaults, case_overrides.get(case_dir.name, {}))
        try:
            print(f"\n>>> 开始处理: {case_dir.name}")
            result = process_one_case(case_dir, case_cfg)
            print(f"    完成: ITC = {result['ITC_W_m2K']:.6e} W/m^2·K")
        except Exception as e:
            result = {
                'case': case_dir.name,
                'status': 'failed',
                'error': str(e),
            }
            print(f"    失败: {e}")
        all_results.append(result)

    summary_csv = root_dir / summary_csv_name
    pd.DataFrame(all_results).to_csv(summary_csv, index=False)
    print(f"\n>>> 全部处理完成，汇总结果已保存到: {summary_csv}")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()
