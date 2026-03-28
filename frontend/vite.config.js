import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

const backendOrigin = process.env.VITE_BACKEND_ORIGIN || 'http://localhost:8000'
const backendWsOrigin =
  process.env.VITE_BACKEND_WS_ORIGIN ||
  backendOrigin.replace(/^http/, 'ws')

export default defineConfig({
  plugins: [react()],
  server: {
    port: 3000,
    proxy: {
      '/api': backendOrigin,
      '/ws': { target: backendWsOrigin, ws: true }
    }
  }
})
