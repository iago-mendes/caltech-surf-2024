\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{letter}[Letter]
\LoadClass[a4paper, 12pt]{article}

% LIBRARIES
	\RequirePackage[top=2cm, bottom=2.3cm, left=2cm, right=2cm, footskip=1.3cm]{geometry}
	\RequirePackage[utf8]{inputenc}
	\RequirePackage{graphicx, caption}
	\RequirePackage{float}
	\RequirePackage{amsmath, amsfonts, amssymb, esint, amsthm, dsfont}
	\RequirePackage{titlesec}
	\RequirePackage{color}
	\RequirePackage{hyperref}
	\RequirePackage{pgf,tikz,pgfplots, tikz-3dplot}
		\usetikzlibrary{arrows}
		\pgfplotsset{compat=1.15}
	\RequirePackage{mathrsfs}
	\RequirePackage{ifthen}
	\RequirePackage{enumitem}
	\RequirePackage{multicol}
	\RequirePackage{subcaption}
	\RequirePackage{fancyhdr}
	\RequirePackage[sorting=none]{biblatex}
	\RequirePackage{etoolbox}


% COLORS
	\definecolor{red}{rgb}{1,0,0}
	\definecolor{green}{rgb}{0,1,0}
	\definecolor{blue}{rgb}{0,0,1}
	\definecolor{link}{HTML}{387494}

	\hypersetup {
		colorlinks=true,
		linkcolor=blue,
		filecolor=magenta,      
		urlcolor=link,
		citecolor=link,
	}

	\newcommand{\red}[1]{\textcolor{red}{#1}}
	\newcommand{\green}[1]{\textcolor{green}{#1}}
	\newcommand{\blue}[1]{\textcolor{blue}{#1}}
	\newcommand{\gray}[1]{\textcolor{gray}{#1}}


% % HEADERS & FOOTERS
	\fancypagestyle{default} {
		\fancyhead{} % clear all header fields
		\renewcommand{\headrule}{\hrulefill}

		\fancyfoot{} % clear all footer fields
		\renewcommand{\footrule}{\hrulefill}
		\fancyfoot[R]{\thepage}
	}
	\pagestyle{default}

	% Show footnotes in footer
	\newcommand{\makefootnotelist}[1] {%
		\parbox{0.8\textwidth} {%
			\footnotesize{%
				\renewcommand*{\do}[1]{##1\\}%
				\dolistcsloop{#1}%
			}%
		}%
	}%
	\newcommand{\fnote}[1] {%
		\footnotemark{}%
		\def\listname{footlist\thepage}%
		\def\n{$^{\the\numexpr\value{footnote}}$}%
		\ifcsdef{\listname}%
			{\listcseadd{\listname}{\n\ #1}}%
			{\csedef{\listname}{}\listcseadd{\listname}{\n\ #1}}%
		\fancypagestyle{fancyfootnote}{%
			\fancyfoot[L]{\makefootnotelist{\listname}}%
			\fancyfoot[R]{\thepage}%
			\fancyfoot[C]{}%
		}%
		\thispagestyle{fancyfootnote}%
	}


% SECTIONS
	\titleformat
		{\section} % command
		{\bf\large} % shape
		{} % label
		{0cm} % separation
		{} % before code
		[] % after code


% COMMANDS

	\renewcommand{\title}[2][a] {
		\begin{center}
			\vspace*{0cm}

			% Title
			{\Large\bf #2}

			\vspace*{0.5cm}
			
			% Subtitle
			\def\argument{#1} \def\default{a}
			\ifx\argument\default
				{}
			\else
				{\large #1}
			\fi

			\vspace*{0.1cm}

			\large
			
			% Author
			Iago Mendes
			
			\vspace*{0.1cm}

			% Date
			\today

			\vspace*{0.75cm}
		\end{center}
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
