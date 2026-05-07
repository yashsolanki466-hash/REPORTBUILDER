export const METAGENOMICS_TEMPLATE = `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Unigenome RNA-seq Report | {{PROJECT_ID}}</title>
    <!-- Fonts & Icons -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    
    <!-- Libraries -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/simple-datatables@9.0.3/dist/style.min.css">
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-zoom@2.2.0"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-zoom@2.2.0"></script>
    <script src="https://unpkg.com/lucide@latest"></script>
    <script src="https://cdn.jsdelivr.net/npm/simple-datatables@9.0.3/dist/umd/simple-datatables.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/apexcharts"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/html2pdf.js/0.10.1/html2pdf.bundle.min.js"></script>

    <style>
        :root {
            /* Unigenome Brand Colors - Official */
            --primary: #2B4570;        /* Navy Blue - Unigenome Brand */
            --primary-light: #4A6FA5;
            --primary-dark: #1A2A45;
            --accent-orange: #FF8C42;  /* Orange - Unigenome Accent */
            --secondary: #10B981;      /* Emerald Green */
            --secondary-light: #6EE7B7;
            --secondary-dark: #059669;
            --accent: #F1F5F9;         /* Light Slate BG - Slightly darker for better contrast */
            --orange: #FF8C42;         /* Unigenome Orange Accent */
            
            /* UI Colors - Enhanced */
            --text-main: #0F172A;
            --text-light: #64748B;
            --text-muted: #94A3B8;
            --border: #E2E8F0;
            --border-light: #F1F5F9;
            --bg-body: #F8FAFC;
            --bg-sidebar: #FFFFFF;
            --bg-card: #FFFFFF;
            --bg-gradient-start: #FAFBFC;
            --bg-gradient-end: #F8FAFC;
            
            /* Status Colors */
            --success: #10B981;
            --success-bg: #D1FAE5;
            --danger: #EF4444;
            --danger-bg: #FEE2E2;
            --warning: #F59E0B;
            --warning-bg: #FEF3C7;

            /* Enhanced Properties */
            --shadow-sm: 0 1px 3px 0 rgb(43 69 112 / 0.08), 0 1px 2px -1px rgb(43 69 112 / 0.06);
            --shadow-md: 0 4px 6px -1px rgb(43 69 112 / 0.1), 0 2px 4px -2px rgb(43 69 112 / 0.08);
            --shadow-lg: 0 10px 15px -3px rgb(43 69 112 / 0.12), 0 4px 6px -4px rgb(43 69 112 / 0.1);
            --shadow-xl: 0 20px 25px -5px rgb(43 69 112 / 0.15), 0 8px 10px -6px rgb(43 69 112 / 0.1);
            --radius-sm: 6px;
            --radius-md: 10px;
            --radius-lg: 16px;
            --radius-xl: 20px;
            --transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            --transition-fast: all 0.15s cubic-bezier(0.4, 0, 0.2, 1);
            --sidebar-width: 280px;
            --sidebar-collapsed-width: 72px;
        }

        .dark {
            --primary: #60A5FA;
            --primary-light: #93C5FD;
            --primary-dark: #2563EB;
            --secondary: #FDBA74;
            --accent: #121212;
            --text-main: #FFFFFF;
            --text-light: #E2E8F0;
            --text-muted: #94A3B8;
            --border: #27272a;
            --border-light: #18181b;
            --bg-body: #000000;
            --bg-sidebar: #050505;
            --bg-card: #0a0a0a;
            --bg-gradient-start: #0F172A; /* Dark blue-ish tint for gradients */
            --success-bg: rgba(16, 185, 129, 0.15);
            --danger-bg: rgba(239, 68, 68, 0.15);
            --warning-bg: rgba(245, 158, 11, 0.15);
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        body {
            font-family: 'Inter', sans-serif;
            background-color: var(--bg-body);
            color: var(--text-main);
            display: flex;
            height: 100vh;
            overflow: hidden;
            transition: background-color 0.3s, color 0.3s;
        }

        /* Modern Scrollbar */
        ::-webkit-scrollbar { width: 6px; height: 6px; }
        ::-webkit-scrollbar-track { background: transparent; }
        ::-webkit-scrollbar-thumb { background: var(--border); border-radius: 10px; }
        ::-webkit-scrollbar-thumb:hover { background: var(--text-muted); }

        /* === SIDEBAR === */
        .sidebar {
            width: var(--sidebar-width);
            background: var(--bg-sidebar);
            border-right: 1px solid var(--border);
            display: flex;
            flex-direction: column;
            flex-shrink: 0;
            z-index: 50;
            transition: width var(--transition);
            overflow: hidden;
            backdrop-filter: blur(10px);
        }
        
        body.sidebar-collapsed .sidebar { width: var(--sidebar-collapsed-width); }
        body.sidebar-collapsed .brand-text { opacity: 0; width: 0; pointer-events: none; }
        body.sidebar-collapsed .nav-group-label { display: none; }
        body.sidebar-collapsed .nav-text { opacity: 0; width: 0; display: none; }
        body.sidebar-collapsed .nav-item { padding: 12px; justify-content: center; }
        body.sidebar-collapsed .nav-item i { margin-right: 0; }
        body.sidebar-collapsed .brand { justify-content: center; padding: 20px 12px; }

        .brand {
            padding: 24px;
            height: 90px;
            display: flex;
            align-items: center;
            border-bottom: 1px solid var(--border);
            /* Removed fixed gradient background for better dark mode compatibility */
            /* background: linear-gradient(135deg, #2B4570 0%, #1A2A45 100%); */
            background: var(--bg-sidebar); 
            white-space: nowrap;
            overflow: hidden;
            transition: var(--transition);
            position: relative;
        }
        
        .brand::before {
            content: '';
            position: absolute;
            inset: 0;
            background: linear-gradient(135deg, transparent 0%, rgba(255, 140, 66, 0.1) 100%);
            pointer-events: none;
        }
        
        .brand::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 0;
            right: 0;
            height: 3px;
            background: linear-gradient(90deg, var(--accent-orange) 0%, var(--primary) 50%, var(--accent-orange) 100%);
            opacity: 0.8;
        }

        .brand-logo { display: flex; align-items: center; gap: 12px; }

        .logo-badge {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            padding: 10px 14px;
            border-radius: 16px;
            background: rgba(255, 255, 255, 0.92);
            border: 1px solid rgba(255, 255, 255, 0.25);
            box-shadow: 0 10px 24px rgba(0, 0, 0, 0.12);
            backdrop-filter: blur(8px);
        }

        .logo-img {
            height: 36px;
            max-width: 190px;
            width: auto;
            object-fit: contain;
            display: block;
        }

        .brand-icon {
            width: 48px;
            height: 48px;
            background: rgba(255, 140, 66, 0.15);
            border: 2px solid rgba(255, 140, 66, 0.3);
            border-radius: 12px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 22px;
            font-weight: 800;
            flex-shrink: 0;
            color: var(--accent-orange);
            backdrop-filter: blur(8px);
            box-shadow: 0 4px 12px rgba(255, 140, 66, 0.2);
        }

        .brand-text { transition: opacity 0.2s; }
        .brand-text h2 { color: white; font-size: 19px; font-weight: 800; letter-spacing: -0.5px; text-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .brand-text p { color: rgba(255,255,255,0.85); font-size: 11px; font-family: 'Inter', sans-serif; font-weight: 500; letter-spacing: 0.3px; }

        .nav-scroll { flex: 1; overflow-y: auto; overflow-x: hidden; padding: 16px 12px; }
        .nav-group-label {
            padding: 16px 12px 8px;
            font-size: 11px;
            text-transform: uppercase;
            letter-spacing: 1px;
            color: var(--text-muted);
            font-weight: 700;
            white-space: nowrap;
        }

        .nav-item {
            display: flex;
            align-items: center;
            padding: 11px 16px;
            margin-bottom: 6px;
            color: var(--text-light);
            text-decoration: none;
            font-size: 14px;
            font-weight: 500;
            cursor: pointer;
            border-radius: var(--radius-md);
            transition: var(--transition-fast);
            white-space: nowrap;
            position: relative;
            overflow: hidden;
        }

        .nav-item::before {
            content: '';
            position: absolute;
            left: 0;
            top: 0;
            height: 100%;
            width: 3px;
            background: linear-gradient(180deg, var(--accent-orange), var(--primary));
            transform: scaleY(0);
            transition: transform 0.2s ease;
            border-radius: 0 2px 2px 0;
        }

        .nav-item:hover { 
            background: linear-gradient(90deg, rgba(255, 140, 66, 0.08) 0%, rgba(43, 69, 112, 0.05) 100%);
            color: var(--primary);
            transform: translateX(4px);
            box-shadow: var(--shadow-sm);
        }
        
        .nav-item.active { 
            background: linear-gradient(90deg, rgba(255, 140, 66, 0.12) 0%, rgba(43, 69, 112, 0.08) 100%);
            color: var(--primary); 
            font-weight: 600;
            box-shadow: var(--shadow-sm);
        }
        
        .nav-item.active::before {
            transform: scaleY(1);
        }
        
        .dark .nav-item.active { background: rgba(59, 130, 246, 0.15); }
        .nav-item i { width: 20px; height: 20px; margin-right: 12px; flex-shrink: 0; stroke-width: 2px; }
        
        /* Progress indicator badge */
        .nav-item .progress-badge {
            margin-left: auto;
            padding: 2px 6px;
            background: var(--success-bg);
            color: var(--success);
            font-size: 10px;
            font-weight: 700;
            border-radius: 10px;
            opacity: 0;
            transform: scale(0.8);
            transition: all 0.2s ease;
        }
        
        .nav-item.completed .progress-badge {
            opacity: 1;
            transform: scale(1);
        }
        
        /* Mobile Menu */
        .mobile-menu-toggle {
            display: none;
            position: fixed;
            bottom: 20px;
            right: 20px;
            width: 56px;
            height: 56px;
            background: var(--primary);
            border-radius: 50%;
            box-shadow: 0 4px 12px rgba(43, 69, 112, 0.4);
            z-index: 100;
            cursor: pointer;
            align-items: center;
            justify-content: center;
            color: white;
            transition: all 0.3s ease;
        }
        
        .mobile-menu-toggle:hover {
            transform: scale(1.1);
            box-shadow: 0 6px 20px rgba(37, 99, 235, 0.6);
        }
        
        /* === RESPONSIVE DESIGN === */
        
        /* Tablet (768px - 1024px) */
        @media (max-width: 1024px) and (min-width: 769px) {
            .sidebar {
                width: 220px;
            }
            .top-header {
                padding: 0 20px;
            }
            .container {
                padding: 0 20px;
            }
            .grid-4 {
                grid-template-columns: repeat(2, 1fr);
            }
        }
        
        /* Mobile (up to 768px) */
        @media (max-width: 768px) {
            .mobile-menu-toggle { 
                display: flex; 
            }
            .sidebar {
                position: fixed;
                left: -100%;
                top: 0;
                height: 100vh;
                width: 280px;
                z-index: 999;
                transition: left 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                box-shadow: 2px 0 12px rgba(0, 0, 0, 0.15);
            }
            .sidebar.mobile-open {
                left: 0;
            }
            .overlay {
                display: none;
                position: fixed;
                inset: 0;
                background: rgba(0, 0, 0, 0.5);
                z-index: 998;
                backdrop-filter: blur(4px);
            }
            .overlay.active {
                display: block;
            }
            .main-content {
                margin-left: 0 !important;
            }
            .top-header {
                padding: 0 16px;
                height: 64px;
            }
            .scroll-area {
                padding: 16px;
            }
            .container {
                padding: 0;
            }
            .header-actions {
                gap: 8px;
            }
            .header-actions .badge {
                display: none;
            }
            .btn {
                padding: 8px 12px;
                font-size: 13px;
            }
            .btn span {
                display: none;
            }
            .page-title {
                font-size: 22px;
            }
            .page-subtitle {
                font-size: 13px;
            }
        }
        
        /* Mobile Portrait (9:16 ratio - optimized for PDF) */
        @media (max-width: 480px) {
            .sidebar {
                width: 100%;
                max-width: 320px;
            }
            .top-header {
                padding: 0 12px;
                height: 56px;
            }
            .scroll-area {
                padding: 12px;
            }
            .stat-box {
                padding: 16px;
            }
            .stat-value {
                font-size: 28px;
            }
            .stat-label {
                font-size: 11px;
            }
            .card-body {
                padding: 16px;
            }
            .page-title {
                font-size: 20px;
            }
            .page-subtitle {
                font-size: 12px;
            }
            table {
                font-size: 11px;
            }
            th, td {
                padding: 10px 12px;
            }
            .mobile-menu-toggle {
                bottom: 16px;
                right: 16px;
                width: 52px;
                height: 52px;
            }
        }

        /* === MAIN CONTENT === */
        .main-content {
            flex: 1;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            position: relative;
            transition: margin-left var(--transition);
        }

        .top-header {
            height: 72px;
            background: rgba(255, 255, 255, 0.95);
            border-bottom: 1px solid var(--border);
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 0 32px;
            flex-shrink: 0;
            backdrop-filter: blur(12px);
            z-index: 40;
            box-shadow: 0 1px 3px rgba(43, 69, 112, 0.05);
        }

        .breadcrumb { display: flex; align-items: center; gap: 10px; font-size: 14px; color: var(--text-muted); }
        .breadcrumb-active { color: var(--text-main); font-weight: 600; }

        .header-actions { display: flex; gap: 12px; align-items: center; }
        
        .btn {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 10px 18px;
            border-radius: var(--radius-md);
            font-size: 13px;
            font-weight: 600;
            cursor: pointer;
            transition: var(--transition-fast);
            border: 1px solid var(--border);
            background: var(--bg-card);
            color: var(--text-main);
            box-shadow: var(--shadow-sm);
        }
        
        .btn:hover { 
            background: var(--accent); 
            border-color: var(--accent-orange); 
            color: var(--primary); 
            transform: translateY(-2px);
            box-shadow: var(--shadow-md);
        }
        .btn:active { transform: translateY(0); }
        
        .btn-primary { 
            background: linear-gradient(135deg, var(--accent-orange) 0%, #FF7A2E 100%);
            color: white; 
            border-color: var(--accent-orange); 
            box-shadow: 0 4px 12px rgba(255, 140, 66, 0.3);
        }
        .btn-primary:hover { 
            background: linear-gradient(135deg, #FF7A2E 0%, var(--accent-orange) 100%);
            color: white; 
            border-color: #FF7A2E;
            box-shadow: 0 6px 16px rgba(255, 140, 66, 0.4);
        }
        .btn-icon { padding: 8px; width: 38px; height: 38px; justify-content: center; border-radius: 50%; }

        .scroll-area { 
            flex: 1; 
            overflow-y: auto; 
            padding: 32px; 
            scroll-behavior: smooth; 
        }
        .container { max-width: 1600px; margin: 0 auto; width: 100%; }

        /* === SECTIONS & CARDS === */
        .section { display: none; animation: slideUp 0.4s cubic-bezier(0.16, 1, 0.3, 1); margin-bottom: 40px; }
        .section.active { display: block; }
        @keyframes slideUp { from { opacity: 0; transform: translateY(20px); } to { opacity: 1; transform: translateY(0); } }

        .page-header { margin-bottom: 32px; }
        .page-title { 
            font-size: 28px; 
            font-weight: 800; 
            background: linear-gradient(135deg, var(--primary) 0%, var(--accent-orange) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 8px; 
            letter-spacing: -0.5px; 
        }
        .page-subtitle { color: var(--text-light); font-size: 15px; max-width: 600px; line-height: 1.5; }

        .card {
            background: var(--bg-card);
            border: 1px solid var(--border);
            border-radius: var(--radius-lg);
            box-shadow: var(--shadow-md);
            margin-bottom: 24px;
            overflow: hidden;
            transition: var(--transition);
            position: relative;
        }
        
        .card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 3px;
            background: linear-gradient(90deg, var(--accent-orange), var(--primary), var(--accent-orange));
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .card:hover::before {
            opacity: 1;
        }
        
        .card:hover {
            box-shadow: var(--shadow-lg);
            transform: translateY(-2px);
        }
        
        .card-header {
            padding: 20px 24px;
            border-bottom: 1px solid var(--border-light);
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: linear-gradient(to bottom, var(--bg-card) 0%, var(--bg-gradient-start) 100%);
        }
        .card-title { 
            font-size: 16px; 
            font-weight: 700; 
            color: var(--primary); 
            display: flex; 
            align-items: center; 
            gap: 10px;
        }
        .card-title i {
            color: var(--accent-orange);
        }
        .card-body { padding: 24px; transition: all 0.3s ease; }
        .card-toggle { background: none; border: none; cursor: pointer; color: var(--text-main); display: flex; align-items: center; justify-content: center; transition: transform 0.3s; }
        .card-toggle:hover { color: var(--primary); }

        .grid-2 { display: grid; grid-template-columns: repeat(2, 1fr); gap: 24px; }
        .grid-4 { display: grid; grid-template-columns: repeat(4, 1fr); gap: 24px; margin-bottom: 32px; }
        @media (max-width: 1280px) { .grid-4 { grid-template-columns: repeat(2, 1fr); } }
        @media (max-width: 768px) { .grid-2, .grid-4 { grid-template-columns: 1fr; } }

        /* === ENHANCED STAT BOXES === */
        .stat-box {
            background: linear-gradient(135deg, var(--bg-card) 0%, var(--bg-gradient-start) 100%);
            border: 1px solid var(--border);
            border-radius: var(--radius-lg);
            padding: 24px;
            position: relative;
            overflow: hidden;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: var(--shadow-md);
        }
        
        .stat-box::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            background: linear-gradient(90deg, var(--accent-orange), var(--primary));
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .stat-box:hover::before {
            opacity: 1;
        }
        
        .stat-box:hover { 
            transform: translateY(-6px) scale(1.02); 
            box-shadow: var(--shadow-xl);
            border-color: var(--accent-orange);
        }
        
        .stat-box::after {
            content: ''; 
            position: absolute; 
            top: -50%; 
            right: -50%; 
            width: 200px; 
            height: 200px;
            background: radial-gradient(circle, var(--warning-bg) 0%, transparent 70%); /* Use variable for accent tint */
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .stat-box:hover::after {
            opacity: 1;
        }
        
        .stat-label { 
            font-size: 13px; 
            text-transform: uppercase; 
            font-weight: 700; 
            color: var(--text-muted); 
            margin-bottom: 12px; 
            letter-spacing: 0.5px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        .stat-value { 
            font-size: 36px; 
            font-weight: 800; 
            background: linear-gradient(135deg, var(--primary) 0%, var(--accent-orange) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -1.5px; 
            line-height: 1;
            position: relative;
            z-index: 1;
        }
        .stat-desc { 
            font-size: 13px; 
            color: var(--text-light); 
            margin-top: 12px; 
            display: flex; 
            align-items: center; 
            gap: 6px; 
            font-weight: 500;
            position: relative;
            z-index: 1;
        }
        
        /* Sparkline mini chart */
        .stat-sparkline {
            width: 60px;
            height: 24px;
            opacity: 0.6;
        }
        
        /* Trend indicator */
        .stat-trend {
            display: inline-flex;
            align-items: center;
            gap: 4px;
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: 700;
        }
        .stat-trend.up {
            background: rgba(16, 185, 129, 0.1);
            color: #059669;
        }
        .stat-trend.down {
            background: rgba(239, 68, 68, 0.1);
            color: #DC2626;
        }
        .stat-trend i {
            width: 12px;
            height: 12px;
        }
        
        /* Animated counter */
        @keyframes countUp {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        .stat-value.animated {
            animation: countUp 0.6s cubic-bezier(0.16, 1, 0.3, 1);
        }

        /* === ENHANCED TABLES === */
        .table-container { 
            overflow-x: auto; 
            border: 1px solid var(--border); 
            border-radius: var(--radius-md); 
            box-shadow: var(--shadow-sm);
            position: relative;
        }
        table { width: 100%; border-collapse: collapse; font-size: 13px; }
        th { 
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
            padding: 14px 20px; 
            text-align: left; 
            font-weight: 700; 
            color: white;
            border-bottom: 2px solid var(--accent-orange);
            white-space: nowrap;
            text-transform: uppercase;
            font-size: 11px;
            letter-spacing: 0.8px;
            position: sticky;
            top: 0;
            z-index: 10;
            cursor: pointer;
            user-select: none;
            transition: all 0.2s ease;
        }
        th:hover {
            background: linear-gradient(135deg, var(--primary-light) 0%, var(--primary) 100%);
        }
        th.sortable::after {
            content: '⇅';
            margin-left: 8px;
            opacity: 0.3;
            font-size: 10px;
        }
        th.sorted-asc::after {
            content: '↑';
            opacity: 1;
            color: var(--primary);
        }
        th.sorted-desc::after {
            content: '↓';
            opacity: 1;
            color: var(--primary);
        }
        td { 
            padding: 14px 20px; 
            border-bottom: 1px solid var(--border-light); 
            color: var(--text-main); 
            transition: all 0.2s ease;
        }
        tr:last-child td { border-bottom: none; }
        tbody tr {
            transition: all 0.2s ease;
        }
        tbody tr:hover {
            background: linear-gradient(90deg, rgba(255, 140, 66, 0.04), rgba(43, 69, 112, 0.03));
            transform: scale(1.005);
            box-shadow: 0 2px 8px rgba(43, 69, 112, 0.08);
        }
        tbody tr:hover td {
            color: var(--primary);
        }
        .text-right { text-align: right; }
        
        /* Cell number formatting */
        .cell-num {
            font-family: 'JetBrains Mono', monospace;
            font-weight: 600;
        }
        
        /* Row selection */
        tbody tr.selected {
            background: rgba(37, 99, 235, 0.1);
            border-left: 3px solid var(--primary);
        }
        
        /* Table actions */
        .table-actions {
            display: flex;
            gap: 8px;
            margin-bottom: 16px;
            align-items: center;
        }
        .table-search {
            flex: 1;
            max-width: 300px;
            padding: 8px 12px;
            border: 1px solid var(--border);
            border-radius: var(--radius-md);
            font-size: 13px;
            transition: all 0.2s ease;
        }
        .table-search:focus {
            outline: none;
            border-color: var(--primary);
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        }
        
        /* === CHART ENHANCEMENTS === */
        .chart-wrapper {
            position: relative;
            margin-bottom: 24px;
        }
        .chart-controls {
            position: absolute;
            top: 16px;
            right: 16px;
            display: flex;
            gap: 8px;
            z-index: 10;
        }
        .chart-btn {
            padding: 8px 14px;
            background: linear-gradient(135deg, white 0%, var(--bg-gradient-start) 100%);
            border: 1px solid var(--border);
            border-radius: var(--radius-md);
            font-size: 12px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            display: flex;
            align-items: center;
            gap: 6px;
            box-shadow: var(--shadow-sm);
            position: relative;
            overflow: hidden;
        }
        
        .chart-btn::before {
            content: '';
            position: absolute;
            inset: 0;
            background: linear-gradient(135deg, rgba(255, 140, 66, 0.1) 0%, rgba(43, 69, 112, 0.05) 100%);
            opacity: 0;
            transition: opacity 0.2s ease;
        }
        
        .chart-btn:hover::before {
            opacity: 1;
        }
        
        .chart-btn:hover {
            background: linear-gradient(135deg, var(--accent) 0%, rgba(255, 140, 66, 0.08) 100%);
            border-color: var(--accent-orange);
            color: var(--primary);
            transform: translateY(-2px);
            box-shadow: var(--shadow-md);
        }
        
        .chart-btn:active {
            transform: translateY(0);
            box-shadow: var(--shadow-sm);
        }
        
        .chart-btn i {
            width: 14px;
            height: 14px;
            position: relative;
            z-index: 1;
        }
        
        .chart-btn span {
            position: relative;
            z-index: 1;
        }
        
        /* Download button specific styling */
        .chart-btn.download-btn {
            background: linear-gradient(135deg, var(--accent-orange) 0%, #FF7A2E 100%);
            color: white;
            border-color: var(--accent-orange);
            box-shadow: 0 2px 8px rgba(255, 140, 66, 0.25);
        }
        
        .chart-btn.download-btn:hover {
            background: linear-gradient(135deg, #FF7A2E 0%, var(--accent-orange) 100%);
            color: white;
            border-color: #FF7A2E;
            box-shadow: 0 4px 12px rgba(255, 140, 66, 0.35);
            transform: translateY(-3px);
        }
        
        .chart-btn.download-btn::before {
            background: linear-gradient(135deg, rgba(255, 255, 255, 0.2) 0%, transparent 100%);
            opacity: 1;
        }
        
        /* Chart legend enhancements */
        .chart-legend {
            display: flex;
            flex-wrap: wrap;
            gap: 16px;
            margin-top: 16px;
            padding: 16px;
            background: var(--accent);
            border-radius: var(--radius-md);
        }
        .legend-item {
            display: flex;
            align-items: center;
            gap: 8px;
            cursor: pointer;
            padding: 6px 12px;
            border-radius: var(--radius-md);
            transition: all 0.2s ease;
        }
        .legend-item:hover {
            background: white;
            box-shadow: var(--shadow-sm);
        }
        .legend-item.disabled {
            opacity: 0.4;
        }
        .legend-color {
            width: 16px;
            height: 16px;
            border-radius: 4px;
            border: 2px solid white;
            box-shadow: 0 0 0 1px var(--border);
        }
        .legend-label {
            font-size: 13px;
            font-weight: 600;
            color: var(--text-main);
        }
        
        /* Crosshair tooltip */
        .chart-crosshair {
            position: absolute;
            pointer-events: none;
            background: rgba(15, 23, 42, 0.95);
            color: white;
            padding: 12px;
            border-radius: var(--radius-md);
            font-size: 12px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
            z-index: 100;
        }
        .chart-crosshair-title {
            font-weight: 700;
            margin-bottom: 6px;
            font-size: 13px;
        }
        .chart-crosshair-value {
            display: flex;
            align-items: center;
            gap: 8px;
            margin-top: 4px;
        }

        /* === MERMAID === */
        .mermaid { 
            width: 100%; 
            display: flex; 
            justify-content: center; 
            min-height: 400px;
        }
        /* Custom Drawflow Node Styling */
        .drawflow .drawflow-node {
            background: transparent;
            border: none;
            box-shadow: none;
            width: auto;
            min-width: 160px;
            padding: 0;
            display: block;
        }
        .drawflow .drawflow-node.selected .drawflow-node-content { 
            border-color: #333; 
            box-shadow: 0 0 0 2px rgba(0,0,0,0.1); 
        }
        
        /* === HEATMAP ENHANCEMENTS === */
        .heatmap-container {
            position: relative;
            background: white;
            border-radius: var(--radius-lg);
            padding: 24px;
            box-shadow: var(--shadow-sm);
        }
        .heatmap-controls {
            display: flex;
            gap: 16px;
            margin-bottom: 24px;
            flex-wrap: wrap;
            align-items: center;
        }
        .heatmap-search {
            flex: 1;
            min-width: 250px;
            padding: 10px 16px;
            border: 1px solid var(--border);
            border-radius: var(--radius-md);
            font-size: 14px;
            transition: all 0.2s ease;
        }
        .heatmap-search:focus {
            outline: none;
            border-color: var(--primary);
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        }
        .heatmap-btn-group {
            display: flex;
            gap: 8px;
        }
        .heatmap-btn {
            padding: 10px 16px;
            background: white;
            border: 1px solid var(--border);
            border-radius: var(--radius-md);
            font-size: 13px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s ease;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        .heatmap-btn:hover {
            background: var(--accent);
            border-color: var(--primary);
            color: var(--primary);
        }
        .heatmap-btn.active {
            background: var(--primary);
            color: white;
            border-color: var(--primary);
        }
        
        /* Dendrogram styling */
        .dendrogram {
            opacity: 0.7;
            transition: opacity 0.2s ease;
        }
        .dendrogram:hover {
            opacity: 1;
        }
        
        /* Gene highlight */
        .gene-highlight {
            outline: 3px solid var(--primary);
            outline-offset: 2px;
            animation: pulse 1s ease-in-out infinite;
        }
        @keyframes pulse {
            0%, 100% { outline-color: var(--primary); }
            50% { outline-color: var(--primary-light); }
        }
        
        /* Clustering info panel */
        .cluster-info {
            position: absolute;
            top: 24px;
            left: 24px;
            background: white;
            padding: 16px;
            border-radius: var(--radius-md);
            box-shadow: var(--shadow-lg);
            border: 1px solid var(--border);
            max-width: 300px;
            z-index: 20;
        }
        .cluster-info-title {
            font-weight: 700;
            font-size: 14px;
            margin-bottom: 8px;
            color: var(--text-main);
        }
        .cluster-info-item {
            font-size: 12px;
            color: var(--text-light);
            margin-top: 4px;
            display: flex;
            justify-content: space-between;
        }

        /* === ENHANCED TABS === */
        .tabs {
            display: flex;
            gap: 4px;
            border-bottom: 2px solid var(--border);
            margin-bottom: 24px;
            overflow-x: auto;
            position: relative;
        }
        .tab-btn {
            padding: 13px 22px;
            background: none;
            border: none;
            color: var(--text-light);
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: var(--transition-fast);
            position: relative;
            white-space: nowrap;
            display: flex;
            align-items: center;
            gap: 8px;
            border-radius: var(--radius-md) var(--radius-md) 0 0;
        }
        .tab-btn::after {
            content: '';
            position: absolute;
            bottom: -2px;
            left: 0;
            right: 0;
            height: 3px;
            background: linear-gradient(90deg, var(--accent-orange), var(--primary));
            transform: scaleX(0);
            transition: transform 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            border-radius: 2px 2px 0 0;
        }
        .tab-btn:hover {
            color: var(--primary);
            background: linear-gradient(to bottom, rgba(255, 140, 66, 0.05), transparent);
        }
        .tab-btn.active {
            color: var(--primary);
            background: linear-gradient(to bottom, rgba(255, 140, 66, 0.08), transparent);
        }
        .tab-btn.active::after {
            transform: scaleX(1);
        }
        .tab-content {
            display: none;
            animation: fadeIn 0.3s ease;
        }
        .tab-content.active {
            display: block;
        }
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        /* Tab badge */
        .tab-badge {
            padding: 2px 8px;
            border-radius: 10px;
            font-size: 10px;
            font-weight: 700;
            background: var(--primary);
            color: white;
        }
        
        /* Keyboard shortcut hint */
        .tab-shortcut {
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 10px;
            font-family: 'JetBrains Mono', monospace;
            background: var(--border-light);
            color: var(--text-muted);
            border: 1px solid var(--border);
        }

        /* === BADGES === */
        .badge { 
            padding: 5px 12px; 
            border-radius: 20px; 
            font-size: 11px; 
            font-weight: 700; 
            text-transform: uppercase; 
            display: inline-flex; 
            align-items: center; 
            gap: 6px;
            box-shadow: var(--shadow-sm);
        }
        .badge-success { background: var(--success-bg); color: #065F46; border: 1px solid rgba(16, 185, 129, 0.2); }
        .badge-danger { background: var(--danger-bg); color: #991B1B; border: 1px solid rgba(239, 68, 68, 0.2); }
        .badge-warning { background: var(--warning-bg); color: #92400E; border: 1px solid rgba(245, 158, 11, 0.2); }

        /* === SECTION NUMBERING (Manual) === */
        /* Auto-numbering removed due to conflict with tab display logic */

        /* === SIGNIFICANCE COLOR CODING FOR TABLES === */
        .sig-high {
            background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%) !important;
            color: #991b1b !important;
            font-weight: 700;
            padding: 4px 8px;
            border-radius: 4px;
            border-left: 3px solid #dc2626;
        }
        
        .sig-medium {
            background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%) !important;
            color: #92400e !important;
            font-weight: 600;
            padding: 4px 8px;
            border-radius: 4px;
            border-left: 3px solid #f59e0b;
        }
        
        .sig-low {
            background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%) !important;
            color: #1e40af !important;
            font-weight: 600;
            padding: 4px 8px;
            border-radius: 4px;
            border-left: 3px solid #3b82f6;
        }
        
        .sig-none {
            background: linear-gradient(135deg, #f3f4f6 0%, #e5e7eb 100%) !important;
            color: #6b7280 !important;
            font-weight: 500;
            padding: 4px 8px;
            border-radius: 4px;
        }
        
        /* Significance legend */
        .sig-legend {
            display: flex;
            gap: 16px;
            padding: 12px 16px;
            background: var(--bg-gradient-start);
            border-radius: var(--radius-md);
            margin-bottom: 16px;
            flex-wrap: wrap;
            border: 1px solid var(--border);
        }
        
        .sig-legend-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 11px;
            font-weight: 600;
        }
        
        .sig-legend-box {
            width: 20px;
            height: 20px;
            border-radius: 4px;
            border: 2px solid white;
            box-shadow: 0 0 0 1px var(--border);
        }
        
        /* Fold change color coding */
        .fc-up {
            color: #dc2626;
            font-weight: 700;
        }
        
        .fc-down {
            color: #2563eb;
            font-weight: 700;
        }
        
        .fc-neutral {
            color: var(--text-muted);
            font-weight: 500;
        }

        /* === DELIVERABLES STACK (Collapsible Horizontal) === */
        .folder-stack { display: flex; flex-direction: column; gap: 16px; }
        
        .folder-card {
            background: var(--bg-card);
            border: 1px solid var(--border);
            border-radius: var(--radius-lg);
            overflow: hidden;
            transition: all 0.2s ease;
        }
        
        .folder-card.active {
            box-shadow: var(--shadow-md);
            border-color: var(--primary-light);
        }
        
        .folder-header {
            padding: 16px 24px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            cursor: pointer;
            background: var(--bg-card);
            transition: background 0.2s;
            user-select: none;
        }
        
        .folder-header:hover {
            background: linear-gradient(to right, rgba(255, 140, 66, 0.03), transparent);
        }
        
        .folder-title-group {
            display: flex;
            align-items: center;
            gap: 16px;
        }
        
        .folder-icon {
            width: 42px; height: 42px;
            border-radius: 12px;
            background: rgba(255, 140, 66, 0.1);
            color: var(--accent-orange);
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
            border: 1px solid rgba(255, 140, 66, 0.2);
        }
        
        .folder-info h4 { margin: 0 0 2px 0; font-size: 15px; font-weight: 700; color: var(--text-main); }
        .folder-info span { font-size: 12px; color: var(--text-muted); font-weight: 500; }
        
        .folder-toggle-icon {
            color: var(--text-muted);
            transition: transform 0.3s ease;
        }
        
        .folder-card.active .folder-toggle-icon {
            transform: rotate(180deg);
            color: var(--primary);
        }
        
        .folder-content {
            display: none;
            padding: 24px;
            background: var(--bg-card); /* Use variable instead of #ffffff */
            border-top: 1px solid var(--border-light);
            animation: slideDown 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }
        
        .folder-content.active {
            display: block;
        }
        
        @keyframes slideDown {
            from { opacity: 0; transform: translateY(-10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .file-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 12px;
        }
        
        .file-item {
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 10px 14px;
            border: 1px solid var(--border-light);
            border-radius: var(--radius-md);
            background: #f8fafc;
            transition: all 0.2s;
        }
        
        .file-item:hover {
            background: var(--bg-body); /* Use variable instead of white */
            border-color: var(--primary-light);
            box-shadow: var(--shadow-sm);
            transform: translateY(-1px);
        }
        
        .file-name-group {
            display: flex;
            align-items: center;
            gap: 10px;
            overflow: hidden;
        }
        
        .file-name-text {
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            font-size: 13px;
            color: var(--text-main);
            font-weight: 500;
        }
        
        .file-action-btn {
            background: white;
            border: 1px solid var(--border);
            border-radius: 6px;
            padding: 6px;
            color: var(--primary);
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: all 0.2s;
            flex-shrink: 0;
        }
        
        .file-action-btn:hover {
            background: var(--primary);
            color: white;
            border-color: var(--primary);
        }

        /* === PRINT & PDF STYLES === */
        @media print {
            @page { size: A4; margin: 20mm; }
            html, body {
                overflow: visible !important;
                height: auto !important;
                background: white !important;
                display: block !important;
            }
            .sidebar, .top-header, .btn, .card-toggle, .nav-item, .header-actions, .mobile-menu-toggle, .overlay { display: none !important; }
            .main-content { margin: 0 !important; width: 100% !important; height: auto !important; overflow: visible !important; display: block !important; }
            .scroll-area { overflow: visible !important; padding: 0 !important; height: auto !important; background: white !important; display: block !important; }
            .container { max-width: 100% !important; margin: 0 !important; padding: 0 !important; display: block !important; }
            .section { display: block !important; page-break-before: always !important; margin-top: 50px !important; }
            .section:first-child { page-break-before: avoid !important; margin-top: 0 !important; }
            .tab-content { display: block !important; opacity: 1 !important; visibility: visible !important; margin-bottom: 30px !important; }
            .card { break-inside: avoid !important; border: 1px solid #eee !important; box-shadow: none !important; margin-bottom: 30px !important; }
            .tab-btn, .tab-bar, .tabs { display: none !important; }
            .card-body, .table-container { display: block !important; height: auto !important; overflow: visible !important; }
            h1, h2, h3, h4 { color: black !important; }
            canvas, .apexcharts-canvas { max-width: 100% !important; height: auto !important; }
        }
        
        /* Mobile PDF Optimization (9:16 ratio) */
        @media print and (max-width: 480px) {
            @page { 
                size: 9in 16in; 
                margin: 15mm; 
            }
            body {
                font-size: 10pt !important;
            }
            .page-title {
                font-size: 18pt !important;
            }
            .stat-box {
                padding: 12pt !important;
            }
            table {
                font-size: 9pt !important;
            }
        }

        /* PDF Export Styling overrides */
        html.is-exporting, 
        body.is-exporting {
            overflow: visible !important;
            height: auto !important;
            background: white !important;
        }

        /* === LATEX-STYLE ACADEMIC EXPORT === */
        .is-exporting {
            font-family: "Times New Roman", Times, serif !important;
            font-size: 11pt !important;
            background: white !important;
            color: black !important;
        }

        .is-exporting .sidebar, 
        .is-exporting .top-header, 
        .is-exporting .card-toggle,
        .is-exporting .header-actions,
        .is-exporting .tab-bar,
        .is-exporting .tabs,
        .is-exporting .btn,
        .is-exporting .comp-toggles {
            display: none !important;
        }

        .is-exporting .main-content { 
            margin: 0 !important; 
            padding: 0 !important;
            width: 100% !important; 
            height: auto !important; 
            overflow: visible !important; 
            display: block !important; 
            position: relative !important;
        }

        .is-exporting .scroll-area { 
            overflow: visible !important; 
            height: auto !important; 
            padding: 20mm !important; /* Standard academic margin */
            background: white !important; 
            display: block !important;
        }

        .is-exporting .section { 
            display: block !important; 
            page-break-before: always; 
            margin-bottom: 2cm !important;
            padding-top: 1cm !important;
            border-top: 0.5pt solid #000;
        }

        .is-exporting .page-title {
            font-size: 24pt !important;
            text-align: center;
            margin-bottom: 2cm !important;
            font-weight: bold;
            text-transform: uppercase;
        }

        .is-exporting .card {
            border: none !important;
            box-shadow: none !important;
            margin-bottom: 1.5cm !important;
            break-inside: avoid !important;
        }

        .is-exporting .card-header {
            background: none !important;
            border-bottom: 1pt solid #000 !important;
            padding: 5pt 0 !important;
            margin-bottom: 10pt !important;
        }

        .is-exporting .card-title {
            font-size: 14pt !important;
            font-weight: bold !important;
            color: black !important;
        }

        .is-exporting .card-body {
            padding: 10pt 0 !important;
            line-height: 1.6 !important;
        }

        .is-exporting table {
            border: 1pt solid black !important;
            border-collapse: collapse !important;
            width: 100% !important;
            margin: 10pt 0 !important;
        }

        .is-exporting th {
            background-color: #f0f0f0 !important;
            border-bottom: 1pt solid black !important;
            color: black !important;
            text-transform: none !important;
            font-weight: bold !important;
        }

        .is-exporting td {
            border-bottom: 0.5pt solid #eee !important;
            color: black !important;
        }

        /* Formal Page Numbering & Indexing Simulation */
        .is-exporting .section::before {
            display: block;
            font-size: 10pt;
            color: #666;
            margin-bottom: 2cm;
            content: "UNIGENOME ANALYSIS REPORT | CONFIDENTIAL";
        }
        
        .is-exporting .tab-content {
            display: block !important;
            opacity: 1 !important;
            visibility: visible !important;
            margin-bottom: 20pt !important;
        }

        .is-exporting .container { 
            max-width: 100% !important; 
            display: block !important; 
            width: 100% !important; 
        }

        /* Clean up Charts for Print */
        .is-exporting canvas, .is-exporting .apexcharts-canvas {
            max-width: 100% !important;
            max-height: 300px !important;
        }

        /* === TABS === */
        
.tab-content { display: none; animation: fadeIn 0.3s ease; }
.tab-content.active { display: block; }
@keyframes fadeIn { from { opacity: 0; } to { opacity: 1; } }
        
.comp-toggles { display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 20px; }
.comp-btn {
    padding: 6px 14px; border: 1px solid var(--border); border-radius: 20px; background: var(--bg-card);
    font-size: 12px; font-weight: 600; cursor: pointer; color: var(--text-light); transition: var(--transition);
}
.comp-btn:hover { border-color: var(--primary); color: var(--primary); }
.comp-btn.active { background: var(--primary); color: white; border-color: var(--primary); }
        .comp-toggles { display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 20px; }
        .comp-btn {
            padding: 8px 16px; 
            border: 1px solid var(--border); 
            border-radius: 20px; 
            background: var(--bg-card);
            font-size: 12px; 
            font-weight: 600; 
            cursor: pointer; 
            color: var(--text-light); 
            transition: var(--transition-fast);
            box-shadow: var(--shadow-sm);
        }
        .comp-btn:hover { 
            border-color: var(--accent-orange); 
            color: var(--primary);
            background: linear-gradient(135deg, rgba(255, 140, 66, 0.08), rgba(43, 69, 112, 0.05));
            transform: translateY(-1px);
            box-shadow: var(--shadow-md);
        }
        .comp-btn.active { 
            background: linear-gradient(135deg, var(--accent-orange), #FF7A2E);
            color: white; 
            border-color: var(--accent-orange);
            box-shadow: 0 4px 12px rgba(255, 140, 66, 0.3);
        }

        /* === OVERLAY === */
        .overlay { 
            position: fixed; inset: 0; background: rgba(0,0,0,0.4); z-index: 998; 
            display: none; opacity: 0; transition: opacity 0.3s; backdrop-filter: blur(4px);
            cursor: pointer;
        }
        .overlay.show, .overlay.active { display: block; opacity: 1; }

        /* === LIGHTBOX & MODAL === */
        .modal-overlay {
            position: fixed; inset: 0; background: rgba(0,0,0,0.85); z-index: 1000;
            display: none; align-items: center; justify-content: center;
            opacity: 0; transition: opacity 0.3s ease; backdrop-filter: blur(5px);
        }
        .modal-overlay.active { display: flex; opacity: 1; }
        
        .lightbox-content {
            max-width: 95%; max-height: 95%;
            position: relative;
            box-shadow: 0 20px 50px rgba(0,0,0,0.5);
            border-radius: 8px; overflow: hidden;
            background: white;
            animation: zoomIn 0.3s cubic-bezier(0.175, 0.885, 0.32, 1.275);
        }
        
        .lightbox-image {
            max-width: 100%; max-height: 90vh; display: block;
        }
        
        .modal-card {
            background: var(--bg-card);
            width: 100%; max-width: 600px;
            border-radius: var(--radius-lg);
            border: 1px solid var(--border);
            box-shadow: var(--shadow-xl);
            display: flex; flex-direction: column;
            max-height: 85vh;
            animation: slideUp 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }
        
        .modal-header {
            padding: 20px 24px;
            border-bottom: 1px solid var(--border);
            display: flex; justify-content: space-between; align-items: center;
        }
        
        .modal-title { font-size: 18px; font-weight: 700; color: var(--text-main); display: flex; align-items: center; gap: 10px; }
        
        .modal-body {
            padding: 24px;
            overflow-y: auto;
            color: var(--text-main);
        }
        
        .modal-close {
            background: none; border: none; color: var(--text-muted); cursor: pointer;
            width: 32px; height: 32px; border-radius: 50%;
            display: flex; align-items: center; justify-content: center;
            transition: all 0.2s;
        }
        .modal-close:hover { background: rgba(0,0,0,0.1); color: var(--text-main); }
        
        .glossary-item {
            margin-bottom: 16px; padding-bottom: 16px;
            border-bottom: 1px solid var(--border-light);
        }
        .glossary-term { font-weight: 700; color: var(--primary); font-size: 14px; margin-bottom: 4px; }
        .glossary-def { font-size: 13px; line-height: 1.5; color: var(--text-light); }
        
        @keyframes zoomIn { from { transform: scale(0.9); opacity: 0; } to { transform: scale(1); opacity: 1; } }
        @keyframes slideUp { from { transform: translateY(20px); opacity: 0; } to { transform: translateY(0); opacity: 1; } }
        
        .scale-bar {
            width: 200px;
            height: 12px;
            background: linear-gradient(to right, #2563eb, #ffffff, #ef4444);
            border-radius: 6px;
        }

        /* ApexCharts Customization */
        .apexcharts-canvas {
            margin: 0 auto;
        }
        .apexcharts-tooltip {
            background: #0f172a !important;
            color: #fff !important;
            border: none !important;
            box-shadow: 0 10px 15px -3px rgba(0,0,0,0.1) !important;
        }

        .mobile-toggle { display: none; background: none; border: none; cursor: pointer; color: var(--text-main); }
        .desktop-toggle { background: none; border: none; cursor: pointer; color: var(--text-main); }

        /* Print adjustments */

        /* Class Code Visuals */
        .class-code-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 16px;
            margin-top: 20px;
        }
        .class-code-item {
            display: flex;
            align-items: center;
            gap: 12px;
            padding: 12px;
            background: var(--bg-body);
            border: 1px solid var(--border);
            border-radius: var(--radius-md);
            transition: all 0.2s;
        }
        .class-code-item:hover {
            border-color: var(--primary);
            box-shadow: var(--shadow-sm);
        }
        .class-code-item.highlight {
            border: 2px solid #22c55e;
            background: #f0fdf4;
            transform: scale(1.02);
            box-shadow: 0 4px 12px rgba(34, 197, 94, 0.15);
        }
        .code-symbol {
            width: 32px;
            height: 32px;
            display: flex;
            align-items: center;
            justify-content: center;
            background: var(--accent);
            color: var(--text-main);
            font-weight: 900;
            font-size: 16px;
            border-radius: 8px;
            flex-shrink: 0;
        }
        .highlight .code-symbol {
            background: #22c55e;
            color: white;
        }
        .code-diagram {
            flex-shrink: 0;
            background: white;
            border-radius: 4px;
            border: 1px solid var(--border-light);
        }
        .code-desc {
            font-size: 12px;
            color: var(--text-light);
            line-height: 1.3;
        }
        .highlight .code-desc { font-weight: 600; color: #166534; }

        /* === INTERACTIVE WORKFLOW === */
        .workflow-container {
            width: 100%;
            min-height: 800px;
            position: relative;
            background: linear-gradient(135deg, #fafbfc 0%, #f1f5f9 100%);
            border-radius: var(--radius-lg);
            overflow: hidden;
            padding: 20px;
        }


        .workflow-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            grid-template-rows: repeat(6, 1fr);
            gap: 20px;
            height: 100%;
            min-height: 760px;
        }


        .workflow-node {
            background: white;
            border: 2px solid #e2e8f0;
            border-radius: var(--radius-lg);
            padding: 16px;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            text-align: center;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            overflow: hidden;
        }


        .workflow-node::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            background: linear-gradient(90deg, var(--primary), var(--accent-orange));
            opacity: 0;
            transition: opacity 0.3s ease;
        }


        .workflow-node:hover {
            transform: translateY(-4px) scale(1.02);
            box-shadow: 0 10px 25px rgba(0, 0, 0, 0.15);
            border-color: var(--primary);
        }

        .workflow-node:hover::before {
            opacity: 1;
        }


        .workflow-node.active {
            border-color: var(--primary);
            background: linear-gradient(135deg, rgba(43, 69, 112, 0.05) 0%, rgba(255, 140, 66, 0.05) 100%);
            box-shadow: 0 0 20px rgba(43, 69, 112, 0.2);
        }


        .workflow-node.process {
            background: linear-gradient(135deg, #deebf7 0%, #f0f8ff 100%);
            border-color: #5b9bd5;
        }


        .workflow-node.input {
            background: linear-gradient(135deg, #e2f0d9 0%, #f0fff0 100%);
            border-color: #70ad47;
        }


        .workflow-node.output {
            background: linear-gradient(135deg, #fff2cc 0%, #fffef0 100%);
            border-color: #ffc000;
        }


        .workflow-node.analysis {
            background: linear-gradient(135deg, #fce4ec 0%, #fdf2f8 100%);
            border-color: #e91e63;
        }


        .workflow-node-icon {
            width: 40px;
            height: 40px;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-bottom: 8px;
            font-size: 18px;
            color: white;
            transition: all 0.3s ease;
        }


        .workflow-node.process .workflow-node-icon {
            background: linear-gradient(135deg, #5b9bd5, #4472c4);
        }


        .workflow-node.input .workflow-node-icon {
            background: linear-gradient(135deg, #70ad47, #548235);
        }


        .workflow-node.output .workflow-node-icon {
            background: linear-gradient(135deg, #ffc000, #e6b800);
        }


        .workflow-node.analysis .workflow-node-icon {
            background: linear-gradient(135deg, #e91e63, #ad1457);
        }

        .workflow-node:hover .workflow-node-icon {
            transform: scale(1.1) rotate(5deg);
        }


        .workflow-node-title {
            font-weight: 600;
            font-size: 13px;
            color: var(--text-main);
            margin-bottom: 4px;
            line-height: 1.2;
        }


        .workflow-node-desc {
            font-size: 11px;
            color: var(--text-muted);
            line-height: 1.3;
        }


        .workflow-connections {
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            pointer-events: none;
            z-index: 1;
        }


        .workflow-tooltip {
            position: absolute;
            background: rgba(0, 0, 0, 0.9);
            color: white;
            padding: 12px 16px;
            border-radius: var(--radius-md);
            font-size: 12px;
            max-width: 250px;
            z-index: 1000;
            pointer-events: none;
            opacity: 0;
            transform: translateY(10px);
            transition: all 0.3s ease;
            backdrop-filter: blur(10px);
        }


        .workflow-tooltip.show {
            opacity: 1;
            transform: translateY(0);
        }


        .workflow-legend {
            display: flex;
            gap: 20px;
            margin-bottom: 20px;
            padding: 16px;
            background: rgba(255, 255, 255, 0.8);
            border-radius: var(--radius-md);
            backdrop-filter: blur(10px);
        }


        .workflow-legend-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 12px;
            color: var(--text-main);
        }


        .workflow-legend-color {
            width: 16px;
            height: 16px;
            border-radius: 4px;
            border: 2px solid white;
            box-shadow: 0 0 0 1px var(--border);
        }


        /* Responsive adjustments */
        @media (max-width: 1200px) {
            .workflow-grid {
                grid-template-columns: repeat(2, 1fr);
                grid-template-rows: repeat(8, 1fr);
            }
        }


        @media (max-width: 768px) {
            .workflow-container {
                padding: 16px;
                min-height: 600px;
            }
            
            .workflow-grid {
                grid-template-columns: 1fr;
                grid-template-rows: repeat(12, 1fr);
                gap: 16px;
            }
            
            .workflow-node {
                padding: 12px;
            }
            
            .workflow-node-title {
                font-size: 12px;
            }
            
            .workflow-node-desc {
                font-size: 10px;
            }
            
            .workflow-legend {
                flex-wrap: wrap;
                gap: 12px;
            }
        }
    </style>
</head>
<body>

    <!-- Overlay -->
    <div class="overlay" id="overlay" onclick="toggleSidebarMobile()"></div>

    <!-- Sidebar -->
    <aside class="sidebar" id="sidebar">
        <div class="brand">
            <div class="brand-logo" style="display: flex; align-items: center; gap: 14px; padding: 8px 0;">
                <div class="logo-badge">
                    <img src="{{UNIGENOME_LOGO_SRC}}" alt="Unigenome" class="logo-img" draggable="false" />
                </div>
            </div>
        </div>

        <nav class="nav-scroll">
            <div class="nav-group-label">Project Overview</div>
            <a class="nav-item active" onclick="showSection('overview')"><i data-lucide="layout-dashboard"></i> <span class="nav-text">Dashboard</span></a>
            <a class="nav-item" onclick="showSection('wet-lab')"><i data-lucide="flask-conical"></i> <span class="nav-text">Wet Lab Analysis</span></a>
            <a class="nav-item" onclick="showSection('workflow')"><i data-lucide="git-branch"></i> <span class="nav-text">Workflow & Methods</span></a>

            <div class="nav-group-label">Analysis Data</div>
            <a class="nav-item" onclick="showSection('data')"><i data-lucide="database"></i> <span class="nav-text">Data Statistics & QC</span></a>
            <a class="nav-item" onclick="showSection('bioinformatics')"><i data-lucide="dna"></i> <span class="nav-text">Bioinformatics Analysis</span></a>
            
            <div class="nav-group-label">Results</div>
            <a class="nav-item" onclick="showSection('deliverables')"><i data-lucide="folder-open"></i> <span class="nav-text">Deliverables</span></a>
            
            <div style="margin-top: auto; padding-top: 20px; border-top: 1px solid rgba(255,255,255,0.1);">
                <a class="nav-item" onclick="toggleGlossary()"><i data-lucide="book"></i> <span class="nav-text">Glossary</span></a>
            </div>
        </nav>
    </aside>

    <!-- Main Content -->
    <main class="main-content">
        <header class="top-header">
            <div style="display: flex; align-items: center; gap: 16px;">
                <button class="mobile-toggle" onclick="toggleSidebarMobile()"><i data-lucide="menu"></i></button>
                <button class="desktop-toggle" onclick="toggleSidebarDesktop()"><i data-lucide="menu"></i></button>
                <div class="breadcrumb">
                    <i data-lucide="home" style="width:16px; height:16px;"></i>
                    <span>/</span>
                    <span id="current-section" class="breadcrumb-active">Dashboard</span>
                </div>
            </div>
            
            <div class="header-actions">
                <span class="badge badge-success" style="padding: 6px 12px; font-size: 11px;">
                    <i data-lucide="check-circle" style="width:12px; height:12px;"></i> QC Passed
                </span>
                <button class="btn btn-icon" onclick="toggleDarkMode()" title="Toggle Dark Mode"><i data-lucide="moon"></i></button>
                <button class="btn btn-primary" onclick="exportToPDF()"><i data-lucide="file-down"></i> Export PDF</button>
            </div>
        </header>

        <div class="scroll-area">
            <div class="container">

                <!-- SECTION: DASHBOARD -->
                <section id="overview" class="section active">
                    <div class="page-header">
                        <h1 class="page-title">1. Whole Transcriptome Analysis</h1>
                        <p class="page-subtitle">Project {{PROJECT_ID}} | {{DATE}} | Submitted to: <strong style="color:var(--primary);">{{CLIENT}} ({{INSTITUTE}})</strong></p>
                    </div>

                    <div class="grid-4">
                        <div class="stat-box">
                            <div class="stat-label">Total Samples</div>
                            <div class="stat-value" id="dash-samples">{{TOTAL_SAMPLES}}</div>
                            <div class="stat-desc"><i data-lucide="server" size="14"></i> {{PLATFORM}}</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-label">Total Data</div>
                            <div class="stat-value" id="dash-data">{{TOTAL_DATA_GB}} <small style="font-size:18px">GB</small></div>
                            <div class="stat-desc"><i data-lucide="hard-drive" size="14"></i> {{READ_LENGTH}}</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-label">Avg Mapping Rate</div>
                            <div class="stat-value" style="color: var(--success);" id="dash-mapping">{{MAPPING_RATE}}</div>
                            <div class="stat-desc"><i data-lucide="check-circle" size="14"></i> High Quality Alignment</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-label">Merged Transcripts</div>
                            <div class="stat-value" id="dash-transcripts">{{MERGED_TRANSCRIPTS}}</div>
                            <div class="stat-desc"><i data-lucide="layers" size="14"></i> StringTie Assembly</div>
                        </div>
                    </div>

                    <div class="grid-2">
                        <!-- Project Details -->
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="file-text"></i> Project Details</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <table class="details-table">
                                    <tr><td width="40%">Service Type</td><td class="text-right"><strong>{{SERVICE_TYPE}}</strong></td></tr>
                                    <tr><td>Platform</td><td class="text-right">{{PLATFORM}}</td></tr>
                                    <tr><td>Read Length</td><td class="text-right">{{READ_LENGTH}}</td></tr>
                                    <tr><td>Organism</td><td class="text-right"><em>{{ORGANISM}}</em> ({{GENOME_BUILD}})</td></tr>
                                </table>
                            </div>
                        </div>

                        <!-- Sample Details -->
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="test-tube"></i> Sample Information</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <table class="details-table">
                                    <tr><td width="40%">Sample Type</td><td class="text-right" id="meta-sampleType">{{SAMPLE_TYPE}}</td></tr>
                                    <tr><td>No. of Samples</td><td class="text-right" id="stat-samples">{{TOTAL_SAMPLES}}</td></tr>
                                    <tr><td>Sample Names</td><td class="text-right" id="meta-sampleNames">{{SAMPLE_NAMES}}</td></tr>
                                    <tr><td>Shipping Condition</td><td class="text-right" id="meta-shipping">{{SHIPPING_CONDITION}}</td></tr>
                                </table>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- SECTION: WET LAB -->
                <section id="wet-lab" class="section">
                    <div class="page-header">
                        <h1 class="page-title">2. Wet Lab Analysis</h1>
                        <p class="page-subtitle">Quality control and library preparation methods</p>
                    </div>
                    
                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="microscope"></i> Methodology</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <div class="grid-2">
                                <div>
                                    <h4 style="font-size:15px; font-weight:700; margin-bottom:8px; color:var(--primary);">1. RNA Isolation & QC</h4>
                                    <p style="font-size:13px; color:var(--text-light); line-height:1.6;">
                                        Total RNA was extracted using the QIAGEN RNeasy mini kit. The quantity was measured using Qubit® 4.0 Fluorometer, and RNA integrity (RIN) was assessed via Agilent 4150 TapeStation to ensure high-quality input material (RIN > 7.0).
                                    </p>
                                </div>
                                <div>
                                    <h4 style="font-size:15px; font-weight:700; margin-bottom:8px; color:var(--primary);">2. Library Preparation</h4>
                                    <p style="font-size:13px; color:var(--text-light); line-height:1.6;">
                                        Sequencing libraries were prepared using the NEBNext® Ultra™ RNA Library Prep Kit. Ribosomal RNA was depleted using the KAPA RNA HyperPrep Kit with RiboErase (HMR) to maximize the coverage of coding and non-coding transcripts.
                                    </p>
                                </div>
                            </div>
                        </div>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title">TapeStation Profiles</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <!-- Simulated Gallery using CSS Grid -->
                            <div id="tapestation-gallery" style="display:grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap:16px;">
                                <div style="background:var(--bg-body); border-radius:12px; height:180px; display:flex; flex-direction:column; align-items:center; justify-content:center; border:1px dashed var(--border);">
                                    <i data-lucide="image" style="color:var(--text-muted); width:32px; height:32px;"></i>
                                    <span style="margin-top:8px; font-size:12px; font-weight:600;">Sample Profile [A]</span>
                                </div>
                                <div style="background:var(--bg-body); border-radius:12px; height:180px; display:flex; flex-direction:column; align-items:center; justify-content:center; border:1px dashed var(--border);">
                                    <i data-lucide="image" style="color:var(--text-muted); width:32px; height:32px;"></i>
                                    <span style="margin-top:8px; font-size:12px; font-weight:600;">Sample Profile [B]</span>
                                </div>
                            </div>
                            <div style="margin-top:16px; font-size:13px; color:var(--success);">
                                <strong><i data-lucide="check"></i> Conclusion:</strong> Libraries showed optimal fragment size distribution (avg ~480-540bp) suitable for 2x150 PE sequencing.
                            </div>
                        </div>
                    </div>
                </section>

                <!-- SECTION: WORKFLOW -->
                <section id="workflow" class="section">
                    <div class="page-header">
                        <h1 class="page-title">3. Bioinformatics Workflow</h1>
                        <p class="page-subtitle">Interactive pipeline visualization and software details</p>
                    </div>

                    <div class="grid-2" style="grid-template-columns: 1.8fr 1fr;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="git-branch"></i> Pipeline Diagram</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body" style="padding:20px; overflow:auto; display:flex; justify-content:center;">
                                <!-- Static SVG Workflow Diagram -->
                                <svg width="800" height="1150" viewBox="0 0 800 1150" xmlns="http://www.w3.org/2000/svg" style="font-family: 'Inter', sans-serif; font-size: 12px;">
                                    <defs>
                                        <marker id="arrow" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto" markerUnits="strokeWidth">
                                            <path d="M0,0 L0,6 L9,3 z" fill="#64748b" />
                                        </marker>
                                    </defs>
                                    
                                    <!-- Styles -->
                                    <style>
                                        .node-input { fill: #E2F0D9; stroke: #70AD47; stroke-width: 1px; }
                                        .node-process { fill: #DEEBF7; stroke: #5B9BD5; stroke-width: 1px; }
                                        .node-output { fill: #FFF2CC; stroke: #FFC000; stroke-width: 1px; }
                                        .link { fill: none; stroke: #64748b; stroke-width: 1.5px; marker-end: url(#arrow); }
                                        .text { font-weight: 500; text-anchor: middle; dominant-baseline: middle; fill: #000; }
                                    </style>

                                    <!-- Row 1 -->
                                    <rect x="50" y="50" width="140" height="50" class="node-input" rx="4" />
                                    <text x="120" y="75" class="text">Reference Genome</text>

                                    <rect x="500" y="50" width="180" height="50" class="node-input" rx="4" />
                                    <text x="590" y="75" class="text">All samples RNA-Seq reads</text>

                                    <!-- Row 2 -->
                                    <rect x="280" y="150" width="120" height="50" class="node-process" rx="4" />
                                    <text x="340" y="175" class="text">STAR</text>

                                    <rect x="500" y="150" width="140" height="50" class="node-output" rx="4" />
                                    <text x="570" y="175" class="text">Reads Alignment</text>

                                    <!-- Row 3 -->
                                    <rect x="50" y="250" width="140" height="60" class="node-input" rx="4" />
                                    <text x="120" y="270" class="text">Reference Genome</text>
                                    <text x="120" y="290" class="text">Annotation</text>

                                    <rect x="280" y="250" width="120" height="50" class="node-process" rx="4" />
                                    <text x="340" y="275" class="text">StringTie</text>

                                    <rect x="500" y="250" width="140" height="50" class="node-output" rx="4" />
                                    <text x="570" y="275" class="text">Assembled Tx</text>

                                    <!-- Row 4 -->
                                    <rect x="280" y="350" width="140" height="50" class="node-process" rx="4" />
                                    <text x="350" y="375" class="text">StringTie Merge</text>

                                    <rect x="500" y="350" width="140" height="50" class="node-output" rx="4" />
                                    <text x="570" y="375" class="text">Merged Tx</text>

                                    <!-- Row 5 (Split) -->
                                    <rect x="200" y="450" width="120" height="50" class="node-process" rx="4" />
                                    <text x="260" y="475" class="text">gff-compare</text>

                                    <rect x="20" y="450" width="140" height="60" class="node-output" rx="4" />
                                    <text x="90" y="470" class="text">Transcripts status</text>
                                    <text x="90" y="490" class="text">vs Reference</text>

                                    <rect x="450" y="450" width="120" height="50" class="node-process" rx="4" />
                                    <text x="510" y="475" class="text">StringTie</text>

                                    <!-- Row 6 -->
                                    <rect x="450" y="550" width="160" height="50" class="node-process" rx="4" />
                                    <text x="530" y="575" class="text">Individual Assembly</text>

                                    <!-- Row 7 -->
                                    <rect x="450" y="650" width="120" height="50" class="node-process" rx="4" />
                                    <text x="510" y="675" class="text">prepDE.py</text>

                                    <!-- Row 8 -->
                                    <rect x="450" y="750" width="120" height="50" class="node-process" rx="4" />
                                    <text x="510" y="775" class="text">edgeR</text>

                                    <!-- Row 9 -->
                                    <rect x="430" y="850" width="160" height="50" class="node-output" rx="4" />
                                    <text x="510" y="875" class="text">Sig. DEGs</text>

                                    <!-- Row 10 (Split) -->
                                    <rect x="280" y="950" width="120" height="50" class="node-process" rx="4" />
                                    <text x="340" y="975" class="text">Blast2GO</text>

                                    <rect x="600" y="950" width="140" height="50" class="node-process" rx="4" />
                                    <text x="670" y="975" class="text">KEGG KAAS</text>

                                    <!-- Row 11 -->
                                    <rect x="280" y="1050" width="120" height="50" class="node-output" rx="4" />
                                    <text x="340" y="1075" class="text">GO Analysis</text>

                                    <rect x="600" y="1050" width="140" height="50" class="node-output" rx="4" />
                                    <text x="670" y="1075" class="text">KEGG Pathways</text>


                                    <!-- CONNECTIONS -->
                                    <!-- Ref Genome -> STAR -->
                                    <path d="M120,100 L120,125 L300,125 L300,144" class="link" />
                                    <!-- Clean Reads -> STAR -->
                                    <path d="M590,100 L590,125 L380,125 L380,144" class="link" />
                                    
                                    <!-- STAR -> Reads Align -->
                                    <path d="M400,175 L494,175" class="link" />
                                    <!-- STAR -> StringTie -->
                                    <path d="M340,200 L340,244" class="link" />
                                    
                                    <!-- Ref Annot -> StringTie -->
                                    <path d="M190,280 L274,280" class="link" />
                                    
                                    <!-- StringTie -> Assembled Tx -->
                                    <path d="M400,275 L494,275" class="link" />
                                    <!-- StringTie -> Merge -->
                                    <path d="M340,300 L340,344" class="link" />
                                    
                                    <!-- Ref Annot -> Merge -->
                                    <path d="M120,310 L120,375 L274,375" class="link" />
                                    
                                    <!-- Merge -> Merged Tx -->
                                    <path d="M420,375 L494,375" class="link" />
                                    
                                    <!-- Merge -> gff-compare -->
                                    <path d="M340,400 L340,425 L260,425 L260,444" class="link" />
                                    <!-- gff-compare -> Tx Status -->
                                    <path d="M200,475 L166,475" class="link" />
                                    
                                    <!-- Merge -> StringTie 2 -->
                                    <path d="M360,400 L360,425 L510,425 L510,444" class="link" />
                                    
                                    <!-- StringTie 2 -> Indiv Assembly -->
                                    <path d="M510,500 L510,544" class="link" />
                                    <!-- Indiv Assembly -> prepDE -->
                                    <path d="M510,600 L510,644" class="link" />
                                    <!-- prepDE -> edgeR -->
                                    <path d="M510,700 L510,744" class="link" />
                                    <!-- edgeR -> Sig DEGs -->
                                    <path d="M510,800 L510,844" class="link" />
                                    
                                    <!-- Sig DEGs -> Blast2GO -->
                                    <path d="M510,900 L510,925 L340,925 L340,944" class="link" />
                                    <!-- Blast2GO -> GO Analysis -->
                                    <path d="M340,1000 L340,1044" class="link" />
                                    
                                    <!-- Sig DEGs -> KEGG KAAS -->
                                    <path d="M510,900 L510,925 L670,925 L670,944" class="link" />
                                    <!-- KEGG KAAS -> KEGG Pathways -->
                                    <path d="M670,1000 L670,1044" class="link" />

                                </svg>
                            </div>
                        </div>

                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="code"></i> Software Details</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="table-container" style="border:none; border-radius:0;">
                                <table>
                                    <thead>
                                        <tr>
                                            <th>Software</th>
                                            <th>Version</th>
                                            <th>Application</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr><td><strong>Trim Galore</strong></td><td>0.6.4</td><td>Adapter and low-quality data removal</td></tr>
                                        <tr><td><strong>STAR</strong></td><td>2.7.10</td><td>Mapping of reads to reference genome</td></tr>
                                        <tr><td><strong>Stringtie</strong></td><td>2.2.1</td><td>Reference based assembly</td></tr>
                                        <tr><td><strong>gffcompare</strong></td><td>0.12.6</td><td>Comparing reference genome gff with merged gff</td></tr>
                                        <tr><td><strong>gffread</strong></td><td>0.12.7</td><td>Fetching fasta sequences from gff</td></tr>
                                        <tr><td><strong>Blastx</strong></td><td>2.30.0+</td><td>Similarity search against NCBI NR</td></tr>
                                        <tr><td><strong>Blast2go cli</strong></td><td>1.5</td><td>GO mapping and annotation</td></tr>
                                        <tr><td><strong>edgeR</strong></td><td>3.6.2</td><td>Differential expression analysis</td></tr>
                                        <tr><td><strong>KEGG KAAS</strong></td><td>Webserver</td><td>Pathway analysis</td></tr>
                                        <tr><td><strong>clusterProfiler</strong></td><td>4.12.6</td><td>Over-representation analysis</td></tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- SECTION: DATA & QC -->
                <section id="data" class="section">
                    <div class="page-header">
                        <h1 class="page-title">4. Data Statistics & QC</h1>
                        <p class="page-subtitle">Raw data metrics and reference genome information</p>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="table"></i> Raw Data Statistics</div>
                            <div style="position: relative;">
                                <i data-lucide="search" style="position: absolute; left: 10px; top: 50%; transform: translateY(-50%); width: 14px; color: var(--text-muted);"></i>
                                <input type="text" id="dataSearch" placeholder="Search samples..." style="padding: 6px 12px 6px 32px; border: 1px solid var(--border); border-radius: 20px; font-size: 13px; background: var(--bg-body); color: var(--text-main);">
                            </div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="table-container">
                            <table id="dataTable">
                                <thead>
                                    <tr>
                                        <th>File Name</th>
                                        <th>Format</th>
                                        <th>Total Reads</th>
                                        <th>Total Bases</th>
                                        <th>GC Content</th>
                                        <th>Q30 %</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <!-- Injected via Script -->
                                </tbody>
                            </table>
                        </div>
                    </div>

                    <div class="grid-4">
                        <div class="stat-box">
                            <div class="stat-label">Genome Size</div>
                            <div class="stat-value">2.85 Gb</div>
                            <div class="stat-desc">GRCh38.p14</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-label">Genes</div>
                            <div class="stat-value">38,558</div>
                            <div class="stat-desc">Reference Annotation</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-label">mRNAs</div>
                            <div class="stat-value">85,576</div>
                            <div class="stat-desc">Transcript Isoforms</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-label">Exons</div>
                            <div class="stat-value">1.32 M</div>
                            <div class="stat-desc">Total Exons</div>
                        </div>
                    </div>
                </section>

                <!-- SECTION: BIOINFORMATICS ANALYSIS -->
                <section id="bioinformatics" class="section">
                    <div class="page-header">
                        <h1 class="page-title">5. Bioinformatics Analysis</h1>
                        <p class="page-subtitle">Microbiome diversity, taxonomy, and composition analysis</p>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="info"></i> Analysis Overview (QIIME2)</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <p style="font-size: 13px; line-height: 1.6; color: var(--text-main);">
                                <strong>QIIME 2 (Quantitative Insights Into Microbial Ecology 2)</strong> is a powerful, extensible, and completely open-source bioinformatics platform for microbiome analysis. In this workflow, raw sequencing data was processed to study the microbiome composition and diversity.
                            </p>
                            <br>
                            <h4 style="font-size: 14px; font-weight: 700; margin-bottom: 8px; color: var(--text-main);">Pipeline Steps:</h4>
                            <ul style="font-size: 13px; line-height: 1.6; color: var(--text-main); padding-left: 20px;">
                                <li><strong>Demultiplexing & QC:</strong> Raw reads were demultiplexed and quality filtered to remove low-quality bases and artifacts.</li>
                                <li><strong>Denoising (DADA2/Deblur):</strong> Sequence variants (ASVs) were identified by modeling and correcting Illumina-sequenced amplicon errors.</li>
                                <li><strong>Taxonomic Classification:</strong> ASVs were classified against a reference database (e.g., Silva, Greengenes) to assign taxonomy.</li>
                                <li><strong>Diversity Analysis:</strong> Alpha and Beta diversity metrics were calculated to assess within-sample and between-sample diversity.</li>
                            </ul>
                        </div>
                    </div>

                    <div class="tabs">
                        <button class="tab-btn active" onclick="switchTab(event, 'bioinformatics', 'taxonomy')">Taxonomy Composition</button>
                        <button class="tab-btn" onclick="switchTab(event, 'bioinformatics', 'alpha-div')">Alpha Diversity</button>
                        <button class="tab-btn" onclick="switchTab(event, 'bioinformatics', 'beta-div')">Beta Diversity</button>
                        <button class="tab-btn" onclick="switchTab(event, 'bioinformatics', 'krona')">Krona Plots</button>
                    </div>

                    <!-- TAXONOMY TAB -->
                    <div id="taxonomy" class="tab-content active">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title">Taxonomic Bar Plot</div>
                                <div style="margin-left: auto; display: flex; gap: 8px;">
                                     <button class="chart-btn" onclick="openLightbox('taxaBarChart')" title="Expand"><i data-lucide="maximize-2"></i></button>
                                     <button class="chart-btn download-btn" onclick="downloadPlot('taxaBarChart', 'taxonomy_barplot.png')"><i data-lucide="download"></i> PNG</button>
                                </div>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 16px;">
                                    Relative abundance of bacterial taxa at different classification levels (Phylum, Genus, etc.). Top 10 most abundant taxa are shown.
                                </p>
                                <div style="height: 400px; width: 100%;">
                                    <canvas id="taxaBarChart"></canvas>
                                </div>
                            </div>
                        </div>

                         <div class="card mt-4">
                            <div class="card-header">
                                <div class="card-title">Taxonomy Table (Top Constituents)</div>
                                <div style="position: relative; margin-left: auto;">
                                    <i data-lucide="search" style="position: absolute; left: 10px; top: 50%; transform: translateY(-50%); width: 14px; color: var(--text-muted);"></i>
                                    <input type="text" id="taxaSearch" placeholder="Search taxa..." style="padding: 6px 12px 6px 32px; border: 1px solid var(--border); border-radius: 20px; font-size: 13px; background: var(--bg-body); color: var(--text-main);">
                                </div>
                            </div>
                            <div class="table-container">
                                <table id="taxaTable">
                                    <thead>
                                        <tr>
                                            <th>Taxon Name</th>
                                            <th>Taxonomic Level</th>
                                            <th>Mean Abundance (%)</th>
                                            <th>Prevalence</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <!-- Injected -->
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>

                    <!-- ALPHA DIVERSITY TAB -->
                    <div id="alpha-div" class="tab-content" style="display: none;">
                        <div class="card">
                             <div class="card-header">
                                <div class="card-title">Alpha Diversity Boxplots</div>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; color: var(--text-main); margin-bottom: 16px;">
                                    <strong>Alpha diversity</strong> measures the species richness and evenness within a single sample. Common metrics include Shannon Index, Observed OTUs, and Chao1.
                                </p>
                                <div class="grid-2">
                                    <div class="card" style="box-shadow:none; border:1px solid var(--border-light);">
                                        <div class="card-header" style="padding:10px;">
                                            <div style="font-size:12px; font-weight:600;">Shannon Diversity</div>
                                             <div style="margin-left: auto; display: flex; gap: 6px;">
                                                <button class="chart-btn" onclick="openLightbox('shannonPlot')" title="Expand"><i data-lucide="maximize-2" style="width:12px; height:12px;"></i></button>
                                            </div>
                                        </div>
                                        <div class="card-body" style="height: 300px;">
                                            <canvas id="shannonPlot"></canvas>
                                        </div>
                                    </div>
                                    <div class="card" style="box-shadow:none; border:1px solid var(--border-light);">
                                        <div class="card-header" style="padding:10px;">
                                            <div style="font-size:12px; font-weight:600;">Observed Species (ASVs)</div>
                                             <div style="margin-left: auto; display: flex; gap: 6px;">
                                                <button class="chart-btn" onclick="openLightbox('observedPlot')" title="Expand"><i data-lucide="maximize-2" style="width:12px; height:12px;"></i></button>
                                            </div>
                                        </div>
                                        <div class="card-body" style="height: 300px;">
                                            <canvas id="observedPlot"></canvas>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- BETA DIVERSITY TAB -->
                    <div id="beta-div" class="tab-content" style="display: none;">
                        <div class="card">
                             <div class="card-header">
                                <div class="card-title">Beta Diversity (PCoA)</div>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; color: var(--text-main); margin-bottom: 16px;">
                                    <strong>Beta diversity</strong> compares the differences between microbial communities from different samples. Principal Coordinates Analysis (PCoA) based on Bray-Curtis or Jaccard distance matrices helps visualize these differences.
                                </p>
                                <div style="height: 500px; width: 100%; display: flex; justify-content: center; align-items: center; background: #f8fafc; border-radius: 8px;">
                                    <!-- Placeholder for PCoA Plot -->
                                     <canvas id="pcoaPlot"></canvas>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- KRONA TAB -->
                    <div id="krona" class="tab-content" style="display: none;">
                         <div class="card">
                             <div class="card-header">
                                <div class="card-title">Interactive Krona Chart</div>
                                <button class="chart-btn" onclick="window.open('krona.html', '_blank')" title="Open Fullscreen" style="margin-left: auto;"><i data-lucide="external-link"></i> Open Fullscreen</button>
                            </div>
                            <div class="card-body" style="height: 600px; padding: 0;">
                                <iframe id="kronaFrame" src="about:blank" style="width: 100%; height: 100%; border: none;"></iframe>
                            </div>
                        </div>
                    </div>
                </section>


                <!-- SECTION: REFERENCE & GFF -->
                <section id="reference" class="section">
                    <div class="page-header">
                        <h1 class="page-title">6. Reference Genome & GFF</h1>
                        <p class="page-subtitle">Details of the reference genome and annotation file used</p>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="dna"></i> Reference Genome Information</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <p style="margin-bottom: 20px; color: var(--text-main); line-height: 1.6; font-size: 13px;">
                                Based on information received from client, the reference genome of <em>Homo sapiens</em> (human) (GCF_000001405.40) and its corresponding GFF file was downloaded from NCBI.
                                <br>
                                <a href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/" target="_blank" style="color: var(--primary); text-decoration: underline;">Source Link</a>
                            </p>

                            <h4 style="font-size: 14px; font-weight: 700; margin-bottom: 12px; color: var(--text-main);">Table 2: Reference genome stats</h4>
                            <div class="table-container">
                                <table>
                                    <thead>
                                        <tr>
                                            <th>Description</th>
                                            <th>Stats</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr><td>Total Chromosomes/Scaffolds</td><td>705</td></tr>
                                        <tr><td>Total genome length (bp)</td><td>3,298,430,636</td></tr>
                                        <tr><td>Mean chromosome/scaffold size</td><td>4,678,625</td></tr>
                                        <tr><td>Max chromosome/scaffold size (bp)</td><td>248,956,422</td></tr>
                                        <tr><td>Min chromosome/scaffolds size (bp)</td><td>970</td></tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="file-text"></i> GFF Annotation Details</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <p style="margin-bottom: 20px; color: var(--text-main); line-height: 1.6; font-size: 13px;">
                                The reference GFF file contains genome annotation. It includes information (locus and description) regarding genes and proteins. The fasta sequence of the reference genome along with the corresponding GFF file is provided with the deliverables in the folder named <strong>"02_reference_genome_and_gff"</strong> for client's reference.
                            </p>

                            <h4 style="font-size: 14px; font-weight: 700; margin-bottom: 12px; color: var(--text-main);">Table 3: Detail Information of GFF file</h4>
                            <div class="table-container">
                                <table>
                                    <thead>
                                        <tr>
                                            <th>Genome elements</th>
                                            <th># Sequences</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr><td>mRNA</td><td>144,447</td></tr>
                                        <tr><td>CDS</td><td>1,836,136</td></tr>
                                        <tr><td>exon</td><td>2,301,289</td></tr>
                                        <tr><td>gene</td><td>47,876</td></tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- SECTION: TRANSCRIPT ASSEMBLY -->
                <section id="assembly" class="section">
                    <div class="page-header">
                        <h1 class="page-title">7. Transcript Assembly</h1>
                        <p class="page-subtitle">StringTie assembly and merging process results</p>
                    </div>

                    <div class="tabs">
                        <button class="tab-btn active" onclick="switchTab(event, 'assembly', 'assembly-main')">Assembly Results</button>
                        <button class="tab-btn" onclick="switchTab(event, 'assembly', 'assembly-comparison')">Comparison & Novel Isoforms</button>
                    </div>

                    <div id="assembly-main" class="tab-content active">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="info"></i> 5.6. Transcript Assembly Methodology</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                        <div class="card-body">
                            <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 16px;">
                                StringTie assembles transcripts from RNA-seq reads that have been aligned to the genome, first grouping the reads into distinct gene loci and then assembling each locus into as many isoforms as are needed to explain the data. Following this, StringTie simultaneously assembles and quantify the final transcripts by using network flow algorithm and starting from most highly abundant transcripts. The reference genome GFF annotation files, containing exon structures of "known" genes, are then used to annotate the assembled transcripts and quantify the expression of known genes as well derive clues if a novel transcript has been found in the sample.
                            </p>
                            <p style="font-size: 13px; line-height: 1.6; color: var(--text-main);">
                                After assembling each sample, the full set of assemblies is passed to StringTie’s merge function, which merges together all the gene structures found in any of the samples. This step is required because transcripts in some of the samples might only be partially covered by reads, and as a consequence only partial versions of them will be assembled in the initial StringTie run. The merge step creates a set of transcripts that is consistent across all samples, so that the transcripts can be compared in subsequent steps.
                            </p>
                        </div>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="image"></i> Figure 11: Merging Transcript Assemblies</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body" style="display:flex; flex-direction:column; align-items:center; background:#fcfcfc;">
                            <!-- SVG Diagram for Transcript Merging -->
                            <svg width="700" height="350" viewBox="0 0 700 350" xmlns="http://www.w3.org/2000/svg" style="font-family: inherit;">
                                <style>
                                    .row-label { font-size: 14px; fill: #334155; }
                                    .exon-green { fill: #22c55e; }
                                    .exon-black { fill: #334155; }
                                    .exon-blue { fill: #2563eb; }
                                    .intron-line { stroke: #94a3b8; stroke-width: 1px; }
                                    .intron-line-bold { stroke: #334155; stroke-width: 1.5px; }
                                    .label-isoform { font-size: 12px; font-style: italic; fill: #64748b; }
                                </style>
                                
                                <!-- Sample 1 -->
                                <text x="20" y="45" class="row-label">Sample 1</text>
                                <line x1="120" y1="40" x2="500" y2="40" class="intron-line" />
                                <rect x="120" y="32" width="60" height="16" class="exon-green" rx="2" />
                                <rect x="300" y="32" width="50" height="16" class="exon-green" rx="2" />
                                <rect x="380" y="32" width="20" height="16" class="exon-green" rx="2" />
                                <rect x="440" y="32" width="60" height="16" class="exon-green" rx="2" />

                                <!-- Sample 2 -->
                                <text x="20" y="95" class="row-label">Sample 2</text>
                                <line x1="150" y1="90" x2="400" y2="90" class="intron-line" />
                                <rect x="150" y="82" width="35" height="16" class="exon-green" rx="2" />
                                <rect x="300" y="82" width="40" height="16" class="exon-green" rx="2" />
                                <rect x="380" y="82" width="12" height="16" class="exon-green" rx="2" />

                                <!-- Sample 3 -->
                                <text x="20" y="145" class="row-label">Sample 3</text>
                                <line x1="320" y1="140" x2="650" y2="140" class="intron-line" />
                                <rect x="320" y="132" width="40" height="16" class="exon-green" rx="2" />
                                <rect x="440" y="132" width="80" height="16" class="exon-green" rx="2" />
                                <rect x="600" y="132" width="50" height="16" class="exon-green" rx="2" />

                                <!-- Sample 4 -->
                                <text x="20" y="195" class="row-label">Sample 4</text>
                                <line x1="130" y1="190" x2="620" y2="190" class="intron-line" />
                                <rect x="125" y="182" width="60" height="16" class="exon-green" rx="2" />
                                <rect x="300" y="182" width="50" height="16" class="exon-green" rx="2" />
                                <rect x="440" y="182" width="80" height="16" class="exon-green" rx="2" />
                                <rect x="600" y="182" width="22" height="16" class="exon-green" rx="2" />

                                <!-- Reference -->
                                <text x="20" y="245" class="row-label" style="font-weight:700;">Reference</text>
                                <text x="20" y="265" class="row-label" style="font-weight:700;">annotation</text>
                                <line x1="120" y1="250" x2="650" y2="250" class="intron-line-bold" />
                                <rect x="120" y="242" width="60" height="16" class="exon-black" rx="2" />
                                <rect x="300" y="242" width="50" height="16" class="exon-black" rx="2" />
                                <rect x="380" y="242" width="20" height="16" class="exon-black" rx="2" />
                                <rect x="440" y="242" width="80" height="16" class="exon-black" rx="2" />
                                <rect x="600" y="242" width="50" height="16" class="exon-black" rx="2" />

                                <!-- Merged Assemblies -->
                                <text x="20" y="315" class="row-label" style="font-weight:700;">Merged</text>
                                <text x="20" y="335" class="row-label" style="font-weight:700;">assemblies</text>
                                
                                <!-- (A) -->
                                <line x1="120" y1="310" x2="650" y2="310" style="stroke: #2563eb; stroke-width: 1px;" />
                                <rect x="120" y="302" width="60" height="16" class="exon-blue" rx="2" />
                                <rect x="300" y="302" width="50" height="16" class="exon-blue" rx="2" />
                                <rect x="380" y="302" width="20" height="16" class="exon-blue" rx="2" />
                                <rect x="440" y="302" width="80" height="16" class="exon-blue" rx="2" />
                                <rect x="600" y="302" width="50" height="16" class="exon-blue" rx="2" />
                                <text x="665" y="315" class="label-isoform">(A)</text>
                                
                                <!-- (B) -->
                                <line x1="125" y1="335" x2="650" y2="335" style="stroke: #2563eb; stroke-width: 1px;" />
                                <rect x="125" y="327" width="55" height="16" class="exon-blue" rx="2" />
                                <rect x="300" y="327" width="50" height="16" class="exon-blue" rx="2" />
                                <rect x="440" y="327" width="80" height="16" class="exon-blue" rx="2" />
                                <rect x="600" y="327" width="50" height="16" class="exon-blue" rx="2" />
                                <text x="665" y="340" class="label-isoform">(B)</text>
                            </svg>
                            <p style="font-size: 12px; color: var(--text-muted); margin-top: 16px; text-align: center; max-width: 600px;">
                                Figure 11: Explanation of merging transcript assemblies using StringTie’s merge function. (Ref: Micheal P et al., 2016).
                            </p>
                        </div>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="table"></i> Table 5: Transcript Assembly Statistics</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="table-container">
                            <table>
                                <thead>
                                    <tr>
                                        <th>Sample Name</th>
                                        <th># assembled transcripts</th>
                                        <th>Total assembled (bp)</th>
                                        <th>Mean transcript size (bp)</th>
                                        <th>Max transcript size (bp)</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr><td><strong>merged.fasta</strong></td><td>6,53,846</td><td>9,830,199,883</td><td>1,503</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-NT-125</td><td>4,84,292</td><td>3,552,472,280</td><td>733.5</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-NT-126</td><td>3,09,175</td><td>6,103,864,499</td><td>1,974.20</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-NT-127</td><td>3,04,981</td><td>5,983,561,143</td><td>1,961.90</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-NT-129</td><td>3,21,143</td><td>6,165,359,962</td><td>1,919.80</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-TT-125</td><td>4,16,107</td><td>4,176,511,182</td><td>1,003.70</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-TT-126</td><td>2,55,210</td><td>5,641,163,378</td><td>2,210.40</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-TT-127</td><td>2,82,706</td><td>5,090,917,59</td><td>1,800.80</td><td>155,532</td></tr>
                                    <tr><td>ICGA-BC-TT-129</td><td>2,69,530</td><td>5,779,757,01</td><td>2,144.40</td><td>155,532</td></tr>
                                </tbody>
                            </table>
                        </div>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="file-check"></i> Output & Naming Notes</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 12px;">
                                The output consists of assembled gene/transcript GTF file for all samples and another GTF file resulting from the merge step as described above. These files have been provided along with the deliverables in the folder named <strong>"03_transcript_assembly_gtf"</strong>. In addition, the transcript sequences derived from these individual assembly is provided in folder <strong>"04_transcript_sequences_fasta"</strong>.
                            </p>
                            <div style="background: var(--bg-body); padding: 12px; border-left: 4px solid var(--primary); border-radius: 4px;">
                                <p style="font-size: 13px; font-weight: 500; color: var(--text-main);">
                                    <strong>Note:</strong> “MSTRG" in output files generated from stringtie stands for "Merged StringTie Transcripts," followed by a unique number. This identifier is created by StringTie merged step to label transcript isoforms that were not previously annotated in the reference genome used.
                                </p>
                            </div>
                        </div>
                        </div>
                    </div>

                    <div id="assembly-comparison" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="git-branch"></i> 5.7. Comparison with Reference Transcripts</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 16px;">
                                    For comparison of assembled transcripts with reference transcripts, <strong>Gffcompare</strong> utility was run taking the reference GTF and the string-tie merged GTF file. This produces an output file, which adds to each transcript a "class code" and the name of the transcript from the reference annotation file to check how the predicted transcripts relate to an annotation file.
                                </p>
                                
                                <h4 style="font-size: 14px; font-weight: 700; margin-bottom: 4px; color: var(--text-main);">Figure 12: Gffcompare Class Codes & Visual Descriptions</h4>
                                <p style="font-size: 12px; color: var(--text-muted); margin-bottom: 20px;">Detailed mapping relationship between assembled transcripts (blue) and reference transcripts (black).</p>
                                
                                <div class="class-code-grid">
                                    <!-- = -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">=</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="20" y="10" width="15" height="10" fill="#334155"/>
                                            <rect x="55" y="10" width="20" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="20" y="25" width="15" height="10" fill="#2563eb"/>
                                            <rect x="55" y="25" width="20" height="10" fill="#2563eb"/>
                                        </svg>
                                        <div class="code-desc">Complete, exact match of intron chain</div>
                                    </div>
                                    <!-- c -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">c</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="15" y="10" width="70" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="30" y="25" width="30" height="10" fill="#2563eb"/>
                                        </svg>
                                        <div class="code-desc">Contained in reference (intron compatible)</div>
                                    </div>
                                    <!-- j -->
                                    <div class="class-code-item highlight">
                                        <div class="code-symbol">j</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="20" y="10" width="15" height="10" fill="#334155"/>
                                            <rect x="55" y="10" width="20" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="100" y2="30" stroke="#22c55e" stroke-width="1.5"/>
                                            <rect x="15" y="25" width="25" height="10" fill="#22c55e"/>
                                            <rect x="55" y="25" width="35" height="10" fill="#22c55e"/>
                                        </svg>
                                        <div class="code-desc"><strong>Novel Isoform:</strong> multi-exon with ≥1 junction match</div>
                                    </div>
                                    <!-- k -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">k</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="15" y="25" width="70" height="10" fill="#2563eb"/>
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="30" y="10" width="30" height="10" fill="#334155"/>
                                        </svg>
                                        <div class="code-desc">Containment of reference (reverse containment)</div>
                                    </div>
                                    <!-- m -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">m</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="20" y="10" width="15" height="10" fill="#334155"/>
                                            <rect x="55" y="10" width="20" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="20" y="25" width="55" height="10" fill="#2563eb"/>
                                        </svg>
                                        <div class="code-desc">Retained intron(s), all introns matched</div>
                                    </div>
                                    <!-- o -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">o</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="40" y="10" width="30" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="20" y="25" width="30" height="10" fill="#2563eb"/>
                                        </svg>
                                        <div class="code-desc">Other same strand overlap with reference exons</div>
                                    </div>
                                    <!-- x -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">x</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="40" y="10" width="30" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#ef4444" stroke-width="1"/>
                                            <rect x="20" y="25" width="30" height="10" fill="#ef4444"/>
                                            <path d="M45,22 L55,22 M50,18 L55,22 L50,26" fill="none" stroke="#ef4444" stroke-width="1"/>
                                        </svg>
                                        <div class="code-desc">Exonic overlap on the opposite strand</div>
                                    </div>
                                    <!-- i -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">i</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="90" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="20" y="10" width="15" height="10" fill="#334155"/>
                                            <rect x="65" y="10" width="15" height="10" fill="#334155"/>
                                            <line x1="10" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="40" y="25" width="20" height="10" fill="#2563eb"/>
                                        </svg>
                                        <div class="code-desc">Fully contained within a reference intron</div>
                                    </div>
                                    <!-- u -->
                                    <div class="class-code-item">
                                        <div class="code-symbol">u</div>
                                        <svg class="code-diagram" width="100" height="40" viewBox="0 0 100 40">
                                            <line x1="10" y1="15" x2="45" y2="15" stroke="#334155" stroke-width="1.5"/>
                                            <rect x="20" y="10" width="15" height="10" fill="#334155"/>
                                            <line x1="55" y1="30" x2="90" y2="30" stroke="#2563eb" stroke-width="1"/>
                                            <rect x="65" y="25" width="15" height="10" fill="#2563eb"/>
                                        </svg>
                                        <div class="code-desc">Unknown, intergenic (no overlap)</div>
                                    </div>
                                </div>
                                <p style="font-size: 11px; color: var(--text-muted); margin-top: 20px;">Reference Source: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml</p>
                            </div>
                        </div>

                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="plus-square"></i> 5.7.1. Identification of Novel Isoform Transcripts</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 16px;">
                                    To identify novel isoform transcripts not present in the reference GTF file, gffcompare utility was ran taking the reference GTF and the string-tie merged GTF file.
                                </p>
                                <div style="background: var(--bg-body); padding: 16px; border-left: 4px solid var(--secondary); border-radius: 4px; margin-bottom: 16px;">
                                    <p style="font-size: 13px; line-height: 1.6; color: var(--text-main);">
                                        <strong>Class code "j"</strong> means that the predicted transcript is a potential <strong>novel isoform</strong> that shares at least one splice junction with a reference transcript.
                                    </p>
                                </div>
                                <p style="font-size: 13px; line-height: 1.6; color: var(--text-main);">
                                    A total of <strong><span id="novel-iso-count"></span> novel isoforms</strong> with the class code "j" were extracted and are provided with the deliverables in a file named <code>novel.isoforms.gtf</code> in the folder <code>"03_transcript_assembly_gff"</code>. 
                                </p>
                                <p style="font-size: 12px; color: var(--text-muted); margin-top: 12px; font-style: italic;">
                                    * A brief description of the column names in the GTF file is provided in the same folder in the file named ReadMe.txt.
                                </p>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- SECTION: DGE -->
                <section id="dge" class="section">
                    <div class="page-header">
                        <h1 class="page-title">8. Differential Expression Analysis</h1>
                        <p class="page-subtitle">Abundance estimation and significance testing</p>
                    </div>

                    <div class="tabs">
                        <button class="tab-btn active" onclick="switchTab(event, 'dge', 'dge-info')">Methodology & Comparisons</button>
                        <button class="tab-btn" onclick="switchTab(event, 'dge', 'dge-stats')">DGE Statistics</button>
                        <button class="tab-btn" onclick="switchTab(event, 'dge', 'dge-plots')">Interactive Plots</button>
                        <button class="tab-btn" onclick="switchTab(event, 'dge', 'dge-heatmap')">Heatmap</button>                    
                    </div>

                    <div id="dge-info" class="tab-content active">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="info"></i> 5.8. Differential Expression Analysis Methodology</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 16px;">
                                    Abundances of the transcripts in all samples were estimated using <strong>StringTie</strong> with the help of merged transcripts generated from merge stringtie step. A python program (<strong>prepDE.py</strong>) was used to extract the read count information directly from the files generated by StringTie. For differential expression analysis, sample comparison was made according to grouping information provided by client.
                                </p>
                                
                                <h4 style="font-size: 14px; font-weight: 700; margin: 20px 0 12px; color: var(--text-main);">Group wise sample information</h4>
                                <div class="table-container">
                                    <table style="font-size: 12px;">
                                        <thead>
                                            <tr><th>Comparison</th><th>Description</th></tr>
                                        </thead>
                                        <tbody>
                                            <tr><td><strong>Comparison 1</strong></td><td>TT-125 (ICGA-BC-TT-125) Vs NT-125 (ICGA-BC-NT-125)</td></tr>
                                            <tr><td><strong>Comparison 2</strong></td><td>TT-126 (ICGA-BC-TT-126) Vs NT-126 (ICGA-BC-NT-126)</td></tr>
                                            <tr><td><strong>Comparison 3</strong></td><td>TT-127 (ICGA-BC-TT-127) Vs NT-127 (ICGA-BC-NT-127)</td></tr>
                                            <tr><td><strong>Comparison 4</strong></td><td>TT-129 (ICGA-BC-TT-129) Vs NT-129 (ICGA-BC-NT-129)</td></tr>
                                            <tr><td><strong>Comparison 5</strong></td><td>TT (All) Vs NT (All) Combined Analysis</td></tr>
                                        </tbody>
                                    </table>
                                </div>
                                <p style="font-size: 11px; color: var(--text-muted); margin-top: 4px; font-style: italic;">Table 6: Sample description and comparisons details.</p>

                                <div style="margin-top: 24px; padding: 16px; background: #f8fafc; border-radius: 8px; border: 1px solid var(--border-light);">
                                    <p style="font-size: 13px; color: var(--text-main);">
                                        Differential gene expression was inferred between samples by applying the R package <strong>edgeR</strong>. It is a Bioconductor package based on negative binomial distribution method. The analysis provides tabular result along with the normalized i.e., Counts Per Million mapped reads (CPM) for each involved sample.
                                    </p>
                                </div>

                                <h4 style="font-size: 14px; font-weight: 700; margin: 24px 0 12px; color: var(--text-main);">Table 7: Column description of edgeR output</h4>
                                <div class="table-container">
                                    <table style="font-size: 12px;">
                                        <thead><tr><th>Column Name</th><th>Feature Identifier</th></tr></thead>
                                        <tbody>
                                            <tr><td><strong>logFC</strong></td><td>log2FoldChange: the logarithm (base 2) of the fold change (Test/Control)</td></tr>
                                            <tr><td><strong>logCPM</strong></td><td>Log10 of CPM value</td></tr>
                                            <tr><td><strong>pval</strong></td><td>pvalue for the statistical significance of this change</td></tr>
                                            <tr><td><strong>FDR</strong></td><td>FDR adjusted pvalue (q-value)</td></tr>
                                        </tbody>
                                    </table>
                                </div>

                                <h4 style="font-size: 14px; font-weight: 700; margin: 24px 0 12px; color: var(--text-main);">Table 8: Criteria used to identify regulated genes</h4>
                                <div class="table-container">
                                    <table style="font-size: 12px;">
                                        <thead><tr><th>Condition</th><th>Status</th></tr></thead>
                                        <tbody>
                                            <tr><td>log2FC > 0</td><td>Up regulated</td></tr>
                                            <tr><td>log2FC < 0</td><td>Down regulated</td></tr>
                                            <tr><td>log2FC > 1 and FDR < 0.05</td><td>Significantly up regulated</td></tr>
                                            <tr><td>log2FC < -1 and FDR < 0.05</td><td>Significantly down regulated</td></tr>
                                        </tbody>
                                    </table>
                                </div>

                                <div style="margin-top: 20px; padding: 12px; background: #eff6ff; border-radius: 6px; border-left: 4px solid #3b82f6;">
                                    <p style="font-size: 12px; color: #1d4ed8; font-weight: 500;">
                                        <strong>Note:</strong> We recommend a minimum of triplicates per group in RNA-Seq experiments. However, as there are no replicates in some studies, we performed DEG analysis using edgeR, which is capable of producing reasonable results in single-sample comparisons.
                                    </p>
                                </div>
                            </div>
                        </div>
                    </div>

                    <div id="dge-stats" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title">DGE Statistics Summary</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="table-container">
                                <table id="dgeSummaryTable">
                                    <thead>
                                        <tr><th>Comparison</th><th>Total DEGs</th><th>Down Reg.</th><th>Up Reg.</th><th>Sig. Down <span class="badge badge-danger">FDR <0.05</span></th><th>Sig. Up <span class="badge badge-success">FDR <0.05</span></th><th>Total Significant</th></tr>
                                    </thead>
                                    <tbody></tbody>
                                </table>
                            </div>
                        </div>

                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="bar-chart"></i> Differential Expression Overview</div>
                                <button class="chart-btn download-btn" onclick="downloadPlot('overviewChart', 'dge_overview.png')" style="margin-left: auto; margin-right: 8px;"><i data-lucide="download"></i> PNG</button>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body" style="height: 300px; position: relative;">
                                <canvas id="overviewChart"></canvas>
                            </div>
                        </div>
                    </div>

                    <div id="dge-plots" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header" style="flex-direction: column; align-items: flex-start; gap: 12px;">
                                <div style="display: flex; width: 100%; justify-content: space-between; align-items: center;">
                                    <div class="card-title">Interactive Plots</div>
                                    <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                                </div>
                                <div class="comp-toggles" id="compToggles"></div>
                            </div>
                            <div class="card-body">
                                <p id="comp-desc" style="font-weight:600; font-size:14px; color:var(--primary); margin-bottom:16px;">Select a comparison...</p>
                                <div class="grid-2">
                                    <div class="card" style="box-shadow:none; border:1px solid var(--border-light);">
                                        <div class="card-header" style="padding:10px; display:flex; justify-content:space-between; align-items:center;">
                                            <div style="font-size:12px; font-weight:600;">MA Plot</div>
                                            <div style="display: flex; gap: 6px;">
                                                <button class="chart-btn" onclick="openLightbox('maPlot')" title="Expand" style="padding: 4px 6px; font-size: 10px;"><i data-lucide="maximize-2" style="width:12px; height:12px;"></i></button>
                                                <button class="chart-btn download-btn" onclick="downloadPlot('maPlot', 'ma_plot.png')" style="padding: 4px 6px; font-size: 10px;"><i data-lucide="download" style="width:12px; height:12px;"></i></button>
                                            </div>
                                        </div>
                                        <div class="card-body" style="height: 400px;">
                                            <canvas id="maPlot"></canvas>
                                        </div>
                                    </div>
                                    <div class="card" style="box-shadow:none; border:1px solid var(--border-light);">
                                        <div class="card-header" style="padding:10px; display:flex; justify-content:space-between; align-items:center;">
                                            <div style="font-size:12px; font-weight:600;">Volcano Plot</div>
                                            <div style="display: flex; gap: 6px;">
                                                <button class="chart-btn" onclick="openLightbox('volcanoPlot')" title="Expand" style="padding: 4px 6px; font-size: 10px;"><i data-lucide="maximize-2" style="width:12px; height:12px;"></i></button>
                                                <button class="chart-btn download-btn" onclick="downloadPlot('volcanoPlot', 'volcano_plot.png')" style="padding: 4px 6px; font-size: 10px;"><i data-lucide="download" style="width:12px; height:12px;"></i></button>
                                            </div>
                                        </div>
                                        <div class="card-body" style="height: 400px;">
                                            <canvas id="volcanoPlot"></canvas>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <div id="dge-heatmap" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="layout-grid"></i> Expression Heatmap (Top DEGs)</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; color: var(--text-main); margin-bottom: 24px;">
                                    Interactive heatmap showing expression patterns of the top 30 most significant genes. Colors indicate relative expression: <span style="color:#ef4444;font-weight:700;">Red (Up)</span>, <span style="color:#2563eb;font-weight:700;">Blue (Down)</span>.
                                </p>
                                <div id="heatmap-apex-container" style="min-height: 500px;">
                                    <!-- Rendered by ApexCharts -->
                                </div>
                            </div>
                        </div>
                    </div>

                    <div id="dge-results" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title">Top Differential Expressed Genes</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="table-container" style="max-height: 500px; overflow-y: auto;">
                                <table id="dgeDetailedTable">
                                    <thead>
                                        <tr><th>Gene/Transcript ID</th><th>Log2FC</th><th>LogCPM</th><th>FDR</th><th>Significance</th></tr>
                                    </thead>
                                    <tbody></tbody>
                                </table>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 12px; color: var(--text-muted);">
                                    Showing top 100 most significant genes sorted by FDR. Full results are available in the <code>05_Differential_Expression</code> folder.
                                </p>
                            </div>
                        </div>
                    </div>
                </section>


                <section id="functional" class="section">
                    <div class="page-header">
                        <h1 class="page-title">9. Functional & Enrichment Analysis</h1>
                        <p class="page-subtitle">Interpretation of significant differential expression</p>
                    </div>

                    <div class="card" style="margin-bottom:24px; border: 1px solid var(--border-light); background: #fdfdfd;">
                        <div class="card-header">
                            <div class="card-title">Select Comparison for Analysis</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <div class="comp-toggles" id="funcCompToggles">
                                <!-- Injected -->
                            </div>
                        </div>
                    </div>

                    <div class="tabs">
                        <button class="tab-btn active" onclick="switchTab(event, 'functional', 'go-analysis')">GO Analysis</button>
                        <button class="tab-btn" onclick="switchTab(event, 'functional', 'kegg-analysis')">KEGG Analysis</button>
                        <button class="tab-btn" onclick="switchTab(event, 'functional', 'pathway-stats')">Pathway Statistics</button>
                        <button class="tab-btn" onclick="switchTab(event, 'functional', 'pathway-example')">Example Pathway</button>
                    </div>

                    <!-- GO ANALYSIS TAB -->
                    <div id="go-analysis-tab" class="tab-content active">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="tag"></i> 5.8.4. Gene Ontology (GO) analysis for significant DEG</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13.5px; line-height: 1.65; color: var(--text-main);">
                                    The Gene Ontology project provides controlled vocabularies of defined terms representing gene product properties. These cover three domains: <strong>Cellular Component</strong>, the parts of a cell or its extracellular environment; <strong>Molecular Function</strong>, the elemental activities of a gene product at the molecular level, such as binding or catalysis; and <strong>Biological Process</strong>, operations or sets of molecular events with a defined beginning and end, pertinent to the functioning of integrated living units: cells, tissues, organs, and organisms. GO was assigned to significant differentially expressed transcripts using Blast2go cli. Single gene can be assigned with multiple GO categories and hence multiple GO terms. The GO domain distribution is shown in table below:
                                </p>
                                
                                <div class="table-title" style="margin-top: 25px;">Table 10: GO Distribution of Significantly Differential expression transcripts</div>
                                <div class="table-container">
                                    <table style="font-size: 12px;" class="table-striped">
                                        <thead>
                                            <tr>
                                                <th>Sample Name</th>
                                                <th>Significant DGE</th>
                                                <th># Seq with GO</th>
                                                <th>Biological Process</th>
                                                <th>Cellular Component</th>
                                                <th>Molecular Function</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            <tr><td>Comparison1</td><td>1,261</td><td>504</td><td>445</td><td>463</td><td>450</td></tr>
                                            <tr><td>Comparison2</td><td>5,147</td><td>2,291</td><td>2,041</td><td>2,127</td><td>2,091</td></tr>
                                            <tr><td>Comparison3</td><td>7,601</td><td>3,008</td><td>2,644</td><td>2,774</td><td>2,677</td></tr>
                                            <tr><td>Comparison4</td><td>4,695</td><td>2,137</td><td>1,877</td><td>1,975</td><td>1,929</td></tr>
                                            <tr><td>Comparison5</td><td>188</td><td>70</td><td>62</td><td>67</td><td>62</td></tr>
                                        </tbody>
                                    </table>
                                </div>
                                <p style="font-size: 12px; color: var(--text-muted); margin-top: 12px;">
                                    Assigned Gene Ontology are also provided in deliverables <code>"06_Significant_DGE_GO"</code>.
                                </p>
                            </div>
                        </div>

                        <div class="card mt-4">
                            <div class="card-header">
                                <div class="card-title">Enriched GO Terms Distribution (Selected Comparison)</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="table-container">
                                <table id="goTable">
                                    <thead><tr><th>Category</th><th>Count</th><th>Description</th></tr></thead>
                                    <tbody><!-- Injected --></tbody>
                                </table>
                            </div>
                        </div>
                    </div>

                    <!-- KEGG ANALYSIS TAB -->
                    <div id="kegg-analysis-tab" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="git-branch"></i> 5.8.5. Pathway analysis for significant DGE</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13.5px; line-height: 1.65; color: var(--text-main);">
                                    Ortholog assignment and mapping of Significant differentially expressed transcripts to the biological pathways were performed using <strong>KEGG automatic annotation server (KAAS)</strong>. Significant differentially expressed transcripts were compared against the KEGG database using <strong>BLASTX</strong> with threshold bit-score value of 60 (default). Pathway analysis was performed using all differentially expressed transcripts. The mapped transcripts represented metabolic pathways of major biomolecules such as carbohydrates, lipids, nucleotides, amino acids, glycans, etc. The mapped transcripts also represented the genes involved in metabolism, genetic information processing, and environmental information processing, cellular processes and organismal systems. Detail of the transcripts following in particular pathway along with the pathway images is given in data deliverables entitled <code>"07_Significant_DGE_pathways"</code>.
                                </p>
                            </div>
                        </div>

                        <div class="card mt-4">
                            <div class="card-header">
                                <div class="card-title">KEGG Enrichment Visualization</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body" style="height: 450px;">
                                <canvas id="keggChart"></canvas>
                            </div>
                        </div>
                    </div>

                    <!-- PATHWAY STATISTICS TAB -->
                    <div id="pathway-stats-tab" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="list"></i> Table 11: KEGG pathway statistics for Significant differentially expressed transcripts</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body">
                                <p style="font-size: 13px; color: var(--text-main); margin-bottom: 24px;">
                                    Pathway statistics for the selected comparison are displayed below. Complete statistics and images for all other combinations have been provided as part of data deliverables <code>"07_Significant_DGE_pathways"</code>.
                                </p>
                                <div class="table-container">
                                    <table id="keggStatsTable" class="table-striped">
                                        <thead>
                                            <tr>
                                                <th>#Level 1</th>
                                                <th>Level 2</th>
                                                <th>Counts</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            <!-- Injected dynamically -->
                                        </tbody>
                                    </table>
                                </div>
                            </div>
                        </div>
                    </div>

                    <!-- EXAMPLE PATHWAY TAB -->
                    <div id="pathway-example-tab" class="tab-content" style="display: none;">
                        <div class="card">
                            <div class="card-header">
                                <div class="card-title"><i data-lucide="image"></i> 5.8.6. Representative Pathway Map</div>
                                <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                            </div>
                            <div class="card-body" style="text-align: center; padding: 20px;">
                                <p style="font-size: 13.5px; line-height: 1.6; color: var(--text-main); margin-bottom: 25px; text-align: left;">
                                    The following image illustrates a representative KEGG pathway map (Glycolysis / Gluconeogenesis), demonstrating how significant differentially expressed genes are integrated into biochemical pathway visualizations.
                                </p>
                                <div style="background: white; padding: 15px; border: 1px solid var(--border-light); border-radius: 12px; box-shadow: var(--shadow-md); display: inline-block;">
                                    <img src="data:image/png;base64,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" alt="Representative Pathway Image" style="max-width: 100%; height: auto; border: 1px solid var(--border-color); border-radius: 8px;">
                                </div>
                                <p style="font-size: 11px; color: var(--text-muted); margin-top: 15px;">
                                    Figure 16: Simplified map illustration of a representative KEGG pathway.
                                </p>
                            </div>
                        </div>
                    </div>
                </section>


                <!-- SECTION: ENRICHMENT ANALYSIS -->
                <section id="enrichment" class="section">
                    <div class="page-header">
                        <h1 class="page-title">Over-representation Analysis</h1>
                        <p class="page-subtitle">Functional enrichment using clusterProfiler</p>
                    </div>

                    <div class="card">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="info"></i> 5.8.6. Over-representation (or enrichment) analysis</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 20px;">
                                Over-representation (enrichment) analysis is a statistical approach used to identify biological functions, pathways, or processes that are significantly over-represented among a given list of genes compared to a reference background. This helps in interpreting large-scale transcriptomic data by linking differentially expressed genes (DEGs) to specific molecular mechanisms and biological pathways.
                            </p>
                            <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 20px;">
                                Functional enrichment analysis was performed using the <strong>clusterProfiler R package</strong> to identify over-represented Gene Ontology (GO) terms and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways among the significantly differentially expressed genes (DEGs). The functions <code>enrichGO</code> and <code>enrichKEGG</code> were applied independently to both up-regulated and down-regulated gene sets.
                            </p>
                            <div style="background: #f8fafc; border: 1px solid var(--border-light); border-radius: 8px; padding: 15px; margin-bottom: 20px;">
                                <p style="font-size: 12px; color: var(--text-main); margin: 0;">
                                    <strong>Parameters:</strong> Multiple hypothesis testing correction was performed using the <strong>Benjamini–Hochberg (BH)</strong> method. GO and KEGG terms with <code>pvalueCutoff = 0.05</code> and <code>qvalueCutoff = 0.2</code> were considered statistically significant.
                                </p>
                            </div>
                        </div>
                    </div>

                    <div class="card mt-4">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="bar-chart-3"></i> Enrichment Visualization</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <p style="font-size: 13px; color: var(--text-main); margin-bottom: 24px;">
                                The results of the enrichment analyses were visualized using dot plots and bar plots, highlighting the most significantly enriched GO terms and KEGG pathways.
                            </p>

                            <!-- Primary Tabs: GO vs KEGG -->
                            <div class="tabs">
                                <button class="tab-btn active" onclick="switchEnrichmentCategory(event, 'go')">GO Enrichment</button>
                                <button class="tab-btn" onclick="switchEnrichmentCategory(event, 'kegg')">KEGG Pathways</button>
                            </div>

                            <!-- GO Enrichment Section -->
                            <div id="go-enrichment-section" class="enrichment-category-section" style="display: block;">
                                <!-- Secondary Tabs: Up vs Down -->
                                <div class="tabs" style="margin-top: 16px; border-top: 1px solid var(--border-light); padding-top: 16px;">
                                    <button class="tab-btn active" onclick="switchEnrichmentDirection(event, 'go', 'up')">Up-regulated</button>
                                    <button class="tab-btn" onclick="switchEnrichmentDirection(event, 'go', 'down')">Down-regulated</button>
                                </div>

                                <!-- GO Up-regulated -->
                                <div id="go-up-section" class="enrichment-direction-section" style="display: block;">
                                    <div id="go-up-dot-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="go-up-dot-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">GO enrichment for up-regulated genes (Dot Plot)</h4>
                                            <div style="width: 100%; height: 500px;">
                                                <canvas id="goDotPlotUp"></canvas>
                                            </div>
                                        </div>
                                    </div>

                                    <div id="go-up-bar-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="go-up-bar-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">GO enrichment for up-regulated genes (Bar Plot)</h4>
                                            <div style="width: 100%; height: 600px;">
                                                <canvas id="goBarPlotUp"></canvas>
                                            </div>
                                        </div>
                                    </div>
                                </div>

                                <!-- GO Down-regulated -->
                                <div id="go-down-section" class="enrichment-direction-section" style="display: none;">
                                    <div id="go-down-dot-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="go-down-dot-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">GO enrichment for down-regulated genes (Dot Plot)</h4>
                                            <div style="width: 100%; height: 500px;">
                                                <canvas id="goDotPlotDown"></canvas>
                                            </div>
                                        </div>
                                    </div>

                                    <div id="go-down-bar-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="go-down-bar-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">GO enrichment for down-regulated genes (Bar Plot)</h4>
                                            <div style="width: 100%; height: 600px;">
                                                <canvas id="goBarPlotDown"></canvas>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            <!-- KEGG Pathways Section -->
                            <div id="kegg-enrichment-section" class="enrichment-category-section" style="display: none;">
                                <!-- Secondary Tabs: Up vs Down -->
                                <div class="tabs" style="margin-top: 16px; border-top: 1px solid var(--border-light); padding-top: 16px;">
                                    <button class="tab-btn active" onclick="switchEnrichmentDirection(event, 'kegg', 'up')">Up-regulated</button>
                                    <button class="tab-btn" onclick="switchEnrichmentDirection(event, 'kegg', 'down')">Down-regulated</button>
                                </div>

                                <!-- KEGG Up-regulated -->
                                <div id="kegg-up-section" class="enrichment-direction-section" style="display: block;">
                                    <div id="kegg-up-dot-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="kegg-up-dot-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">KEGG pathway enrichment for up-regulated genes (Dot Plot)</h4>
                                            <div style="width: 100%; height: 500px;">
                                                <canvas id="keggDotPlotUp"></canvas>
                                            </div>
                                        </div>
                                    </div>

                                    <div id="kegg-up-bar-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="kegg-up-bar-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">KEGG pathway enrichment for up-regulated genes (Bar Plot)</h4>
                                            <div style="width: 100%; height: 600px;">
                                                <canvas id="keggBarPlotUp"></canvas>
                                            </div>
                                        </div>
                                    </div>
                                </div>

                                <!-- KEGG Down-regulated -->
                                <div id="kegg-down-section" class="enrichment-direction-section" style="display: none;">
                                    <div id="kegg-down-dot-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="kegg-down-dot-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">KEGG pathway enrichment for down-regulated genes (Dot Plot)</h4>
                                            <div style="width: 100%; height: 500px;">
                                                <canvas id="keggDotPlotDown"></canvas>
                                            </div>
                                        </div>
                                    </div>

                                    <div id="kegg-down-bar-tab" class="plot-type-content" style="display: block;">
                                        <div style="display: flex; flex-direction: column; align-items: center; background: white; padding: 20px; border-radius: 8px; margin-top: 16px;">
                                            <h4 id="kegg-down-bar-title" style="font-size: 13px; font-weight: 700; color: var(--text-main); margin-bottom: 20px;">KEGG pathway enrichment for down-regulated genes (Bar Plot)</h4>
                                            <div style="width: 100%; height: 600px;">
                                                <canvas id="keggBarPlotDown"></canvas>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            <p style="font-size: 11px; color: var(--text-muted); margin-top: 20px; font-style: italic;">
                                Note: Plots for only one comparison have been shown above, however plots for upregulated genes as well as all other combinations have been provided as a part of data deliverables <code>"08_Significant_DGE_Enrichment/"</code>.
                            </p>
                        </div>
                    </div>

                    <div style="margin-top: 24px; padding: 16px; background: #fff1f2; border-radius: 8px; border-left: 4px solid #f43f5e;">
                        <p style="font-size: 13px; color: #9f1239; font-weight: 600;">⚠️ Exception Note:</p>
                        <p style="font-size: 12px; color: #be123c; margin-top: 4px;">
                            No enrichment was obtained for the following combinations: <strong>C5_up_KEGG</strong>.
                        </p>
                    </div>
                </section>

                <!-- SECTION: DELIVERABLES -->
                <section id="deliverables" class="section">
                    <div class="page-header">
                        <h1 class="page-title">10. Deliverables</h1>
                        <p class="page-subtitle">Structure of the output data package</p>
                    </div>

                    <div class="folder-grid" id="folder-grid">
                        <!-- Injected -->
                    </div>

                    <div class="card mt-4" id="uploaded-files-card" style="display: none;">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="file-text"></i> Uploaded Files</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <div class="uploaded-files-grid" id="uploaded-files-grid">
                                <!-- Injected -->
                            </div>
                        </div>
                    </div>

                    <div class="card mt-4" id="ftp-access-card" style="display: none;">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="database"></i> FTP Server Access</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <div id="ftp-login-section">
                                <div class="ftp-login-form" style="max-width: 400px; margin: 0 auto; padding: 20px; border: 1px solid var(--border); border-radius: 8px; background: var(--bg-body);">
                                    <h4 style="text-align: center; margin-bottom: 20px; color: var(--text-main);">FTP Server Login</h4>
                                    <div style="margin-bottom: 16px;">
                                        <label style="display: block; margin-bottom: 4px; font-weight: 600; font-size: 14px; color: var(--text-main);">User ID</label>
                                        <input type="text" id="ftp-username" placeholder="Enter your user ID" style="width: 100%; padding: 8px 12px; border: 1px solid var(--border); border-radius: 4px; font-size: 14px;">
                                    </div>
                                    <div style="margin-bottom: 20px;">
                                        <label style="display: block; margin-bottom: 4px; font-weight: 600; font-size: 14px; color: var(--text-main);">Password</label>
                                        <input type="password" id="ftp-password" placeholder="Enter your password" style="width: 100%; padding: 8px 12px; border: 1px solid var(--border); border-radius: 4px; font-size: 14px;">
                                    </div>
                                    <button id="ftp-login-btn" onclick="connectToFtp()" style="width: 100%; padding: 10px; background: var(--primary); color: white; border: none; border-radius: 6px; font-weight: 600; cursor: pointer; margin-bottom: 16px;">
                                        Connect to FTP Server
                                    </button>
                                    <div id="ftp-status" style="text-align: center; font-size: 12px; color: var(--text-muted);"></div>
                                </div>
                                <div id="ftp-files-section" style="display: none;">
                                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px;">
                                        <h4 style="color: var(--text-main);">Available Files</h4>
                                        <button onclick="disconnectFtp()" style="padding: 6px 12px; background: #ef4444; color: white; border: none; border-radius: 4px; font-size: 12px; cursor: pointer;">Disconnect</button>
                                    </div>
                                    <div id="ftp-files-grid" style="display: grid; gap: 8px;">
                                        <!-- FTP files will be populated here -->
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <div class="card mt-4">
                        <div class="card-header">
                            <div class="card-title"><i data-lucide="folder-tree"></i> Directory Tree Overview</div>
                            <button class="card-toggle" onclick="toggleCard(this)"><i data-lucide="chevron-up"></i></button>
                        </div>
                        <div class="card-body">
                            <pre id="tree-view" style="font-family: 'JetBrains Mono', monospace; font-size: 12px; background: #f8fafc; padding: 20px; border-radius: 8px; border: 1px solid var(--border); overflow-x: auto; color: var(--text-main); margin: 0; line-height: 1.5;"></pre>
                        </div>
                    </div>
                </section>

            </div>
        </div>
    </main>

    <!-- Floating Mobile Menu Button -->
    <button class="mobile-menu-toggle" onclick="toggleSidebarMobile()" aria-label="Toggle Menu">
        <i data-lucide="menu"></i>
    </button>

    <!-- Data Injection Script (At the bottom to not block render) -->
    <!-- Lightbox Overlay -->
    <div class="modal-overlay" id="lightbox-overlay" onclick="closeLightbox()">
        <div class="lightbox-content" onclick="event.stopPropagation()">
            <button class="modal-close" onclick="closeLightbox()" style="position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.8); z-index: 10;"><i data-lucide="x"></i></button>
            <img id="lightbox-image" class="lightbox-image" src="" alt="Fullscreen Chart">
        </div>
    </div>

    <!-- Glossary Modal -->
    <div class="modal-overlay" id="glossary-modal" onclick="toggleGlossary()">
        <div class="modal-card" onclick="event.stopPropagation()">
            <div class="modal-header">
                <div class="modal-title"><i data-lucide="book-open" style="color:var(--primary);"></i> Glossary of Terms</div>
                <button class="modal-close" onclick="toggleGlossary()"><i data-lucide="x"></i></button>
            </div>
            <div class="modal-body" id="glossary-content">
                <!-- Injected via JS -->
            </div>
        </div>
    </div>

    <script>
        // --- 0. UI Functions ---
        window.showSection = function(id) {
            document.querySelectorAll('.section').forEach(el => el.classList.remove('active'));
            document.getElementById(id).classList.add('active');
            
            // Update Active Nav
            document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
            // Simple match
            const navs = document.querySelectorAll('.nav-item');
            for(let nav of navs) {
                if(nav.getAttribute('onclick').includes(id)) {
                    nav.classList.add('active');
                    document.getElementById('current-section').textContent = nav.textContent.trim();
                }
            }
            // Close sidebar on mobile
            if(window.innerWidth < 768) {
                const sb = document.getElementById('sidebar');
                const ov = document.getElementById('overlay');
                if (sb) sb.classList.remove('mobile-open');
                if (ov) ov.classList.remove('active', 'show');
                document.body.style.overflow = '';
            }
        };

        window.toggleSidebarMobile = function() {
            const sb = document.getElementById('sidebar');
            const ov = document.getElementById('overlay');
            if (!sb || !ov) return;
            const isOpen = sb.classList.toggle('mobile-open');
            ov.classList.toggle('active', isOpen);
            ov.classList.toggle('show', isOpen);
            document.body.style.overflow = isOpen ? 'hidden' : '';
        };

        window.toggleSidebarDesktop = function() {
            document.body.classList.toggle('sidebar-collapsed');
        };

        window.toggleDarkMode = function() {
            document.body.classList.toggle('dark');
            // Re-render charts for color adaptation if needed
        };

        window.toggleCard = function(btn) {
            const card = btn.closest('.card');
            // Find content: could be card-body or table-container
            const content = card.querySelector('.card-body') || card.querySelector('.table-container');
            
            if(content) {
                if(content.style.display === 'none') {
                    content.style.display = 'block';
                    btn.style.transform = 'rotate(0deg)';
                } else {
                    content.style.display = 'none';
                    btn.style.transform = 'rotate(180deg)';
                }
            }
        };

        window.toggleFolder = function(header) {
            const card = header.closest('.folder-card');
            const content = card.querySelector('.folder-content');
            
            if(content) {
                if(content.classList.contains('active')) {
                    content.classList.remove('active');
                    card.classList.remove('active');
                } else {
                    content.classList.add('active');
                    card.classList.add('active');
                }
            }
        };

        // --- LIGHTBOX & GLOSSARY FUNCTIONS ---
        window.openLightbox = function(canvasId) {
            const canvas = document.getElementById(canvasId);
            if(!canvas) return;
            
            const overlay = document.getElementById('lightbox-overlay');
            const imgTarget = document.getElementById('lightbox-image');
            
            // Create a temporary white background version for better visibility
            const tempCanvas = document.createElement('canvas');
            tempCanvas.width = canvas.width;
            tempCanvas.height = canvas.height;
            const ctx = tempCanvas.getContext('2d');
            ctx.fillStyle = '#ffffff';
            ctx.fillRect(0, 0, tempCanvas.width, tempCanvas.height);
            ctx.drawImage(canvas, 0, 0);
            
            imgTarget.src = tempCanvas.toDataURL('image/png');
            overlay.classList.add('active');
            document.body.style.overflow = 'hidden';
        };
        
        window.closeLightbox = function() {
            document.getElementById('lightbox-overlay').classList.remove('active');
            document.body.style.overflow = '';
        };
        
        window.toggleGlossary = function() {
            const el = document.getElementById('glossary-modal');
            if(el.classList.contains('active')) {
                el.classList.remove('active');
                document.body.style.overflow = '';
            } else {
                el.classList.add('active');
                document.body.style.overflow = 'hidden';
            }
        };

        // --- GLOSSARY DATA ---
        const glossaryTerms = {
            "FDR (False Discovery Rate)": "A statistical method used in multiple hypothesis testing to correct for random events that falsely appear significant. An FDR of 0.05 implies that 5% of significant results may be false positives.",
            "LogCPM (Log Counts Per Million)": "A normalized measure of gene expression. It represents the log2-transformed count of reads mapped to a gene per million reads sequenced, allowing for comparison across samples with different sequencing depths.",
            "TPM (Transcripts Per Million)": "A normalization method for RNA-seq that accounts for gene length and sequencing depth. TPM is useful for comparing the proportion of transcripts within a sample.",
            "Log2 Fold Change (Log2FC)": "A measure of how much a gene's expression has changed between two conditions. A value of 1 means expression has doubled (up-regulated), while -1 means it has halved (down-regulated).",
            "P-value": "A measure of the probability that an observed difference could have occurred just by random chance. Lower impact values (typically < 0.05) suggest the result is statistically significant.",
            "GO (Gene Ontology)": "A system for classifying gene functions into three categories: Biological Process, Molecular Function, and Cellular Component.",
            "KEGG (Kyoto Encyclopedia of Genes and Genomes)": "A database resource for understanding high-level functions and utilities of the biological system, such as the cell, the organism and the ecosystem, from molecular-level information.",
            "DGE (Differential Gene Expression)": "The process of statistical analysis to discover quantitative changes in expression levels between experimental groups.",
            "PCA (Principal Component Analysis)": "A dimensionality reduction technique used to visualize the variation present in a dataset. In RNA-seq, it helps to see how similar samples are to each other based on their gene expression profiles."
        };
        
        // Populate Glossary
        document.addEventListener('DOMContentLoaded', () => {
            const container = document.getElementById('glossary-content');
            if(container) {
                container.innerHTML = Object.entries(glossaryTerms).map(([term, def]) => \`
                    <div class="glossary-item">
                        <div class="glossary-term">\${term}</div>
                        <div class="glossary-def">\${def}</div>
                    </div>
                \`).join('');
            }
        });

        // --- Download Individual Plot Function (High Quality White Background) ---
        window.downloadPlot = function(canvasId, filename) {
            const originalCanvas = document.getElementById(canvasId);
            if (!originalCanvas) {
                console.error('Canvas not found:', canvasId);
                return;
            }
            
            // Create a temporary canvas to composite the white background
            const tempCanvas = document.createElement('canvas');
            const ctx = tempCanvas.getContext('2d');
            
            // Set dimensions to match original (or scale up for better quality if needed)
            const scale = 2; // Verify if 2x improves quality for the user, usually 1x or 2x is good. Using original w/h for now but handled via direct dimension copy
            tempCanvas.width = originalCanvas.width;
            tempCanvas.height = originalCanvas.height;
            
            // 1. Fill with white background
            ctx.fillStyle = '#ffffff';
            ctx.fillRect(0, 0, tempCanvas.width, tempCanvas.height);
            
            // 2. Draw the original chart
            ctx.drawImage(originalCanvas, 0, 0);
            
            // 3. Add a footer/watermark (Optional premium touch)
            ctx.font = '12px Inter, sans-serif';
            ctx.fillStyle = '#64748b'; // muted blue-grey
            ctx.textAlign = 'right';
            ctx.fillText('Generated by Unigenome', tempCanvas.width - 15, tempCanvas.height - 15);

            // Create download link
            const link = document.createElement('a');
            link.download = filename || (canvasId + '_plot.png');
            link.href = tempCanvas.toDataURL('image/png', 1.0);
            link.click();
            
            // Show success notification
            console.log('Downloaded with white background: ' + link.download);
        };

        // --- Significance Color Coding Helper ---
        window.applySigColorCoding = function(value, type = 'pvalue') {
            // Apply significance classes based on p-value or FDR
            // Returns appropriate CSS class
            const numValue = parseFloat(value);
            
            if (isNaN(numValue)) return 'sig-none';
            
            if (type === 'pvalue' || type === 'fdr') {
                if (numValue < 0.001) return 'sig-high';      // p < 0.001 - Highly significant
                if (numValue < 0.01) return 'sig-medium';     // p < 0.01 - Medium significance
                if (numValue < 0.05) return 'sig-low';        // p < 0.05 - Low significance
                return 'sig-none';                             // Not significant
            }
            
            return 'sig-none';
        };

        // --- Apply Fold Change Color Coding ---
        window.applyFCColorCoding = function(value, threshold = 1) {
            const numValue = parseFloat(value);
            
            if (isNaN(numValue)) return 'fc-neutral';
            
            if (numValue > threshold) return 'fc-up';      // Up-regulated
            if (numValue < -threshold) return 'fc-down';   // Down-regulated
            return 'fc-neutral';                            // Not significant
        };

        // --- Format Significance Value ---
        window.formatSigValue = function(value) {
            const numValue = parseFloat(value);
            if (isNaN(numValue)) return value;
            
            // Scientific notation for very small values
            if (numValue < 0.0001) {
                return numValue.toExponential(2);
            }
            
            // Fixed decimal for normal values
            return numValue.toFixed(4);
        };

        window.switchTab = function(event, sectionId, tabId) {
            const section = document.getElementById(sectionId);
            if(!section) return;
            
            // Hide all tab contents in this section
            const contents = section.querySelectorAll('.tab-content');
            contents.forEach(c => {
                c.style.display = 'none';
                c.classList.remove('active');
            });
            
            // Show selected
            const target = document.getElementById(tabId + '-tab') || document.getElementById(tabId);
            if(target) {
                target.style.display = 'block';
                target.classList.add('active');
            }
            
            // Active button
            const btn = event.currentTarget;
            if(btn) {
                const bar = btn.parentElement;
                bar.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
            }
        };

        // Enrichment Analysis Nested Tab Functions
        window.switchEnrichmentCategory = function(event, category) {
            // Hide all category sections
            document.querySelectorAll('.enrichment-category-section').forEach(sec => {
                sec.style.display = 'none';
            });
            
            // Show selected category
            const targetSection = document.getElementById(category + '-enrichment-section');
            if(targetSection) {
                targetSection.style.display = 'block';
            }
            
            // Update active button
            const btn = event.currentTarget;
            if(btn) {
                const bar = btn.parentElement;
                bar.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
            }
        };

        window.switchEnrichmentDirection = function(event, category, direction) {
            // Hide all direction sections within this category
            const categorySection = document.getElementById(category + '-enrichment-section');
            if(!categorySection) return;
            
            categorySection.querySelectorAll('.enrichment-direction-section').forEach(sec => {
                sec.style.display = 'none';
            });
            
            // Show selected direction
            const targetSection = document.getElementById(category + '-' + direction + '-section');
            if(targetSection) {
                targetSection.style.display = 'block';
            }
            
            // Update active button
            const btn = event.currentTarget;
            if(btn) {
                const bar = btn.parentElement;
                bar.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
            }
        };

        window.exportToPDF = function() {
            const btn = document.querySelector('[onclick="exportToPDF()"]');
            const originalContent = btn.innerHTML;
            
            // 1. Preparation
            btn.innerHTML = '<i data-lucide="loader-2" class="animate-spin"></i> Finalizing PDF...';
            btn.style.pointerEvents = 'none';
            
            // Ensure we are at the top so capture starts correctly
            window.scrollTo(0, 0);
            document.querySelector('.scroll-area').scrollTo(0,0);
            
            document.documentElement.classList.add('is-exporting');
            document.body.classList.add('is-exporting');
            
            const element = document.querySelector('.main-content');
            const projectID = window.REPORT_DATA?.metadata?.projectID || 'Report';
            
            const opt = {
                margin:       [10, 5, 10, 5],
                filename:     'Unigenome_Report_' + projectID + '.pdf',
                image:        { type: 'jpeg', quality: 0.95 },
                html2canvas:  { 
                    scale: 1, // Critical: Scale 1 is required for very long documents (30+ + pages) to avoid canvas limits
                    useCORS: true, 
                    logging: false,
                    letterRendering: true,
                    scrollY: 0
                },
                jsPDF:        { unit: 'mm', format: 'a4', orientation: 'portrait' },
                pagebreak:    { mode: ['css', 'legacy'] }
            };

            // 2. Generation with extra delay for rendering stability
            setTimeout(() => {
                html2pdf().set(opt).from(element).save().then(() => {
                    document.documentElement.classList.remove('is-exporting');
                    document.body.classList.remove('is-exporting');
                    btn.innerHTML = originalContent;
                    btn.style.pointerEvents = 'auto';
                    if(window.lucide) window.lucide.createIcons();
                }).catch(err => {
                    console.error('PDF Export failed:', err);
                    alert('Export failed. Please try again or use the browser print (Ctrl+P) option.');
                    document.documentElement.classList.remove('is-exporting');
                    document.body.classList.remove('is-exporting');
                    btn.innerHTML = originalContent;
                    btn.style.pointerEvents = 'auto';
                });
            }, 2000);
        };

        // --- 1. Init Icons ---
        document.addEventListener('DOMContentLoaded', () => {
            if(window.lucide) window.lucide.createIcons();
        });
    </script>
</body>
</html>`;

