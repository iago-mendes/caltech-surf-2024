\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{letter}[Letter]
\LoadClass[a4paper, 12pt]{article}

\RequirePackage[top=2cm, bottom=2cm, left=2cm, right=2cm]{geometry}
\RequirePackage[utf8]{inputenc}
\RequirePackage{graphicx, caption}
\RequirePackage{float}
\RequirePackage{amsmath, amsfonts, amssymb, esint, amsthm, dsfont}
\RequirePackage{titlesec}
\RequirePackage{color}
\RequirePackage{hyperref}
\RequirePackage{pgf,tikz,pgfplots, tikz-3dplot}
\RequirePackage{mathrsfs}
\RequirePackage{ifthen}
\RequirePackage{enumitem}
\RequirePackage{multicol}
\RequirePackage{subcaption}

\usetikzlibrary{arrows}
\pgfplotsset{compat=1.15}

\hypersetup {
	colorlinks=true,
	linkcolor=blue,
	filecolor=magenta,      
	urlcolor=blue,
	citecolor=blue,
}

% COLORS
\definecolor{red}{rgb}{1,0,0}
\definecolor{green}{rgb}{0,1,0}
\definecolor{blue}{rgb}{0,0,1}
\definecolor{ccqqqq}{rgb}{0.8,0,0}
\newcommand{\red}[1]{\textcolor{red}{#1}}
\newcommand{\green}[1]{\textcolor{green}{#1}}
\newcommand{\blue}[1]{\textcolor{blue}{#1}}
\newcommand{\gray}[1]{\textcolor{gray}{#1}}

\newcommand{\setPageStyle}[1] {
	\newpagestyle{mystyle} {
		% \sethead{}{}{}
		\headrule
		\footrule
			\setfoot{}{\thepage}{}
	}

	\pagestyle{mystyle}
	\thispagestyle{mystyle}
}

\newcommand{\header}[1] {
	\title{\vspace{-20px}#1}
	\author{Iago Mendes}
	\date{\today}
	
	\maketitle

	\setPageStyle{#1}
}

\newcommand{\instructions}[1]{
	\gray{#1} \\
}

\newcommand{\image}[2][0.5] {
	\begin{figure}[H]
		\centering
		\includegraphics[width=#1\paperwidth]{#2}
	\end{figure}
}

% optional parameters: height
% required parameters: image1, label1, image2, label2, description
\newcommand{\twoImages}[6][125px] {
	\begin{figure}[H]
		\centering
		\begin{subfigure}{.4\textwidth}
			\centering
			\includegraphics[height=#1]{#2}
			\caption{}
			\label{#3}
		\end{subfigure}
		\begin{subfigure}{.4\textwidth}
			\centering
			\includegraphics[height=#1]{#4}
			\caption{}
			\label{#5}
		\end{subfigure}
		\caption{#6}
	\end{figure}
}
