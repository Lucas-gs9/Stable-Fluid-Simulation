# Real-time Fluid Simulation

A real-time 2D fluid simulation based on the Stable-Fluids algorithm, featuring RGB color rendering and interactive mouse controls.

## Overview

This project implements Jos Stam's Stable-Fluids method for real-time fluid dynamics simulation. It supports colorful fluid rendering through RGB density fields and allows users to interact with the fluid using mouse gestures.

## Features

- **RGB Color Fluid**: Three-channel density fields for vibrant color mixing
- **Real-time Interaction**: Mouse drag to inject fluid and apply forces
- **Background Velocity Fields**: Horizontal flow, rising smoke, and vortex presets
- **Infinite Domain**: Fluid re-enters from opposite boundaries instead of disappearing

## Build

This project uses [xmake](https://xmake.io) for building.

## Credits

This project is based on the course framework[PKU-VCL 2025](https://gitee.com/pku-vcl/vci-2025 "the VCI-2025 course framework from Peking University Visual Computing Lab"), the Peking University Visual Computing Lab. The core Stable-Fluids implementation was developed independently as part of the coursework.

## Reference
[Jos Stam, "Stable Fluids", SIGGRAPH 1999](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf)

```bash
xmake
xmake run final_project
