// main.js — the FloodA5 globe viewer with the CesiumJS engine.
import { startViewer } from '../../viewer-common/js/app.js';
import { CesiumEngine } from './engine.js';

startViewer(new CesiumEngine());
